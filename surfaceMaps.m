% This code reads in a patient's CT scan and dose grid, 
% puts them in image coordiante space, generates a surface map of
% depth-weighted HU and dose, and finds the surface normal 
% vectors parallel to the camera line

%% READ IN CT AND DOSE FILES AND RESAMPLE 
tic
temp = pwd; %directory of .dcm files
files = dir(temp);
RDfiles = {};
for f = 1:size(files,1)
    if startsWith(files(f).name, 'RD')
        RDfiles = [RDfiles;files(f).name];
    elseif startsWith(files(f).name, 'RP')
        RPfile = files(f).name;
    end
end

% Read CT and header info -------------------

ctlist = dir(fullfile(temp, 'ct')); % moved all ct dicoms into a subdirectory called "CT"
ctslices= cell2struct(struct2cell(ctlist),fieldnames(ctlist));
i = 1;
for c = 1:size(ctlist)
    if startsWith(ctlist(c).name, 'CT') == 1
        ctslices(i) = ctlist(c);
        i = i+1;
    end
end
ctslices = ctslices(1:i);
ctinfo = dicominfo(fullfile(temp, 'ct',ctslices(1).name));  
offset = ctinfo.RescaleIntercept;
ctvol_in = medicalVolume(fullfile(temp, 'ct'));
%ctvol_in = tiffreadVolume('Pt31ctstack.tif'); %some stacks need to be read in with FIJI and saved as a tiff stack. i.e.: test = tiffreadVolume('pt31ctstack.tif');
ctimpospat = ctvol_in.VolumeGeometry.Position(1,:); %image position patient, x,y position of top left pixel 
%ctimpospat = ctinfo.ImagePositionPatient; 
ctres= ctvol_in.VoxelSpacing;
%ctres = [ctinfo.PixelSpacing', ctinfo.SliceThickness] ;
ctvol = double(squeeze(dicomreadVolume(fullfile(temp, 'ct'))));
%ctvol = double(squeeze(ctvol_in));
ctvol = ctvol + offset;
 
% Read dose and header info -------------------
for i = 1:numel(RDfiles)

    doseinfo = dicominfo(fullfile(temp,RDfiles{i}));
    dosescalefactor = doseinfo.DoseGridScaling;
    
    doseimpospat = doseinfo.ImagePositionPatient;
    
    doseres = [doseinfo.PixelSpacing; ctres(3)];
    dosevol = double(squeeze(dicomread(fullfile(temp,RDfiles{i}))));
    dosevol = dosevol*dosescalefactor*100/fx; %dose per fraction
    
    if i == 1
        dosesum = dosevol;
    elseif i > 1
        dosesum = dosesum + dosevol;
    end
end


% Interpolate z direction to 1 mm slice thickness ------------------
ctvol = imresize3(ctvol, [size(ctvol,1), size(ctvol,2), 3*size(ctvol,3)]);

ctres(3) = ctres(3)/3;

%x, y = 2.5mm res, slice thickness same as CT
% 3 is the slice thickness, may need to change to make more general
dosesum = imresize3(dosesum, [size(dosesum,1), size(dosesum,2), 3*size(dosesum,3)]);
doseres(3) = doseres(3)/3;


% Get CT meshgrid ------------------ %CT coord. space
%top left to top right
xx = ctimpospat(1): ctres(1): ctimpospat(1)+ctres(1)*(size(ctvol,2)-1);
yy = ctimpospat(2): ctres(2): ctimpospat(2)+ctres(2)*(size(ctvol,1)-1);
zz = ctimpospat(3): ctres(3): ctimpospat(3)+ctres(3)*(size(ctvol,3)-1);

%make iso 0, 0, 0
xx = xx - iso(1);
yy = yy - iso(2);
zz = zz - iso(3);

[ctXX,ctYY,ctZZ] = meshgrid(xx,yy,zz); 

ctMesh.ctXX = ctXX;
ctMesh.ctYY = ctYY;
ctMesh.ctZZ = ctZZ;

% Get Dose meshgrid ------------------
xx = doseimpospat(1): doseres(1): doseimpospat(1)+doseres(1)*(size(dosesum,2)-1);
yy = doseimpospat(2): doseres(2): doseimpospat(2)+doseres(2)*(size(dosesum,1)-1);
zz = doseimpospat(3): doseres(3): doseimpospat(3)+doseres(3)*(size(dosesum,3)-1);

xx = xx - iso(1);
yy = yy - iso(2);
zz = zz - iso(3);

[dXX,dYY,dZZ] = meshgrid(xx,yy,zz);

dMesh.dXX = dXX;
dMesh.dYY = dYY;
dMesh.dZZ = dZZ;

% Interpolate dose onto ct volume
dose = interp3(dXX,dYY,dZZ,dosesum, ctXX,ctYY,ctZZ, 'cubic',0);

%% GENERATE SURFACE MAPS ORIENTED IN IMAGE SPACE FROM THE CAMERA VIEWING ANGLE

% Crop CT -----------------
def_offset = -1000;
% Set HU values to air around the bench, find slice where couch starts
ctvol_temp = ctvol;
ctvol(400:end,:,:)= def_offset;
%ctvol(:,:,380:end)= def_offset;
%ctvol(:, :, 1:100)= def_offset;
%ctvol(:, :, 1:end) = ctvol(:, :, end:-1:1); %needed if read in as .tif
%from FIJI
%
% Threshold for skin segmentation
HUthres = -150;  

% Define mesh size but maybe use the actual CT coords
[mm,nn,pp] = size(ctvol);

[X,Y,Z] = meshgrid(1:nn, 1:mm, 1:pp); %index coords, interpolating to 1 mm

% Creating a surface from ctvol using HUthresh
surf_iso = isosurface(ctMesh.ctXX,ctMesh.ctYY,ctMesh.ctZZ,ctvol, HUthres); %iso coords
surf_ind = isosurface(X,Y,Z,ctvol, HUthres); %index coords

% collect vertices
v_iso = surf_iso.vertices;
v_ind = surf_ind.vertices;

% Normal vectors pointing outward from the surface
N_ind = isonormals(X,Y,Z,ctvol,v_ind, 'reverse');

% View isosurface
figure;
patch('Faces',surf_iso.faces,'Vertices',surf_iso.vertices,'FaceColor','flat','FaceVertexCData',[0 1 1],'EdgeColor','none');
set(gca,'DataAspectRatio',[1,1,1]);
light 

% Surface weighted dose and HU values -------------------

doseColors = zeros(size(v_ind,1),1); %NX1
ctColors = zeros(size(v_ind,1),1);

f = waitbar(0,'Please wait...');
Nn = numel(ctColors);

d1 = 1; % Begin sampling at this location (in mm), assuming 1 mm spacing
d2 = 5; % End sampling at this location inside the patient

for i=1:Nn
    % Getting the matrix indices of each vertex involved (in isocenter
    % space)
    coords = round(v_ind(i,:));
    xx = coords(1);
    yy = coords(2);
    zz = coords(3);
    % [~, xx, ~] = ind2sub(size(ctMesh.ctXX), find(coords(1) == round(ctMesh.ctXX))); xx = mode(xx);
    % [yy, ~ , ~] = ind2sub(size(ctMesh.ctXX), find(coords(2) == round(ctMesh.ctYY))); yy = mode(yy);
    % [~ , ~ , zz] = ind2sub(size(ctMesh.ctZZ), find(coords(3) == round(ctMesh.ctZZ))); zz = mode(zz);

    vec = -1*N_ind(i,:);
    vec = vec/norm(vec); % Normalized inverse vector to follow for sampling

    sampDose = zeros(numel(d1:d2),1);
    sampCT = zeros(numel(d1:d2),1);

    k = 1;
    for j = d1:d2 % change sample depth, loop through sampling depth
        veci = round(j*vec);
        xxi = xx+veci(1);
        yyi = yy+veci(2);
        zzi = zz+veci(3);

        % To avoid going outside of volume
        if(xxi <= 0 || xxi > max(nn(:))) xxi = 1; end
        if(yyi <= 0 || yyi > max(mm(:)))  yyi = 1; end
        if(zzi <= 0 || zzi > max(pp(:))) zzi = 1; end

        % Query to dose for the specific index
       % sampDose(k) = dose(yyi,xxi,zzi);
         sampCT(k) = ctvol(yyi,xxi,zzi);
         k = k+1;

    end

    % Take the average dose of the following, but leave out NaN values.
    % could do exp. sampling here if wanted
   % doseColors(i)   =   mean(sampDose(:),'omitnan');
    ctColors(i)     =   mean(sampCT(:),'omitnan');

    % Updating the percentages
    if mod(i,5000) == 0
        perc = i/Nn;
        waitbar(perc,f,['Please wait (',num2str(100*perc,2),'%) ...'])
    end
end

close(f)

% Color the faces with the surface dose or HU values-------------------------

%close all;
f1 = figure;
hold on;
p = patch('Faces',surf_iso.faces,'Vertices',surf_iso.vertices);
p.FaceVertexCData = ctColors; %or could be CT
p.FaceColor = 'interp';
p.EdgeColor = 'none';

ax = gca;
set(ax,'DataAspectRatio',[1 1 1])
xlabel('x');
ylabel('y');
zlabel('z');
camup('manual')
clim([-200 200])
colormap(bone)

%Keene left cam
% camup([-0.428322444538922 -0.8614915768571899 0.2727125712326344]);
% campos([1967.1420305273687	-1336.511224453499 -1246.6280925228862]);
% camtarget([13.85196350956221 26.92738968124104 -7.400638824133694]);
% camva(10.19051586894112)

%Moffitt Right Cam
camup([0.4675053886614372 -0.8623565422631594 0.19437053683222144]);
campos([-1615.4530911003208 -1068.4646556278083 -754.6003412378782]);
camtarget([-1.7801874301633234 -17.336446178554752 27.650470285178017]);
camva(11.337405759636438)

camproj('perspective')

cam = plotCamera(AbsolutePose=worldPose,Opacity=0.3,size=50);

% Generate point cloud from mesh vertices --------------------------

%figure('Name','PointCloud','Position',[50 50 1200 1920]);
pointCloud_dose = pointCloud(v_iso);
pointCloud_dose.Intensity = ctColors;
%pcshow(pointCloud_dose,'BackgroundColor',[0 0 0],'ColorSource','Intensity');
% colormap(bone)
% axis off 

% Use external settings to set up camera
% camup([-0.428322444538922 -0.8614915768571899 0.2727125712326344]);
% campos([1967.1420305273687	-1336.511224453499 -1246.6280925228862]);
% camtarget([13.85196350956221 26.92738968124104 -7.400638824133694]);
% camva(10.19051586894112)

%Moffitt Right Cam
% camup([0.4675053886614372 -0.8623565422631594 0.19437053683222144]);
% campos([-1615.4530911003208 -1068.4646556278083 -754.6003412378782]);
% camtarget([-1.7801874301633234 -17.336446178554752 27.650470285178017]);
% camva(11.337405759636438)


%% Remove "invisible" points ----------------------------------------
% see this reference for built in method: https://doi.org/10.1145/1276377.1276407
% PDF: https://www.weizmann.ac.il/math/ronen/sites/math.ronen/files/uploads/katz_tal_basri_-_direct_visibility_of_point_sets.pdf

%need to adjust radius
[pointCloud_dose_vis,visInd] = removeHiddenPoints(pointCloud_dose,...
    [-1615.4530911003208 -1068.4646556278083 -754.6003412378782], RadiusScale=5);

figure('Name','PointCloud','Position',[50 50 1200 1920]);
pcshow(pointCloud_dose_vis, 'BackgroundColor',[0 0 0],'ColorSource','Intensity');
colormap(gray)
axis off 

% Use external settings to set up camera
camup([0.4675053886614372 -0.8623565422631594 0.19437053683222144]);
campos([-1615.4530911003208 -1068.4646556278083 -754.6003412378782]);
camtarget([-1.7801874301633234 -17.336446178554752 27.650470285178017]);
camva(11.337405759636438)


v_iso_2D = world2img(pointCloud_dose_vis.Location,camExtrinsics,camIntrinsics);
%v_iso_2D = world2img(ptCloud2.Location,camExtrinsics,camIntrinsics);
%figure;
% scatter(v_iso_2D(:,1),v_iso_2D(:,2),4,ptCloud2.Intensity,'filled');
% colormap(jet)
% axis equal
% set(gca,'Ydir','reverse')
% xlim([1,1920]);
% ylim([1,1200]);
p2 = pointCloud_dose_vis; p2.Intensity = doseColors(visInd);
maxDose = max(p2.Intensity);


% Interpolate visible point could to image matrix

% F = scatteredInterpolant(v_iso_2D(:,1), v_iso_2D(:,2), pointCloud_dose_vis.Intensity,'linear','none');
% rows = 1200;
% columns = 1920;
% [xGrid, yGrid] = meshgrid(1:columns, 1:rows);
% xq = xGrid(:);
% yq = yGrid(:);
% vq = F(xq, yq);
% ctImage = reshape(vq, rows, columns);
% figure; imagesc(ctImage)
% colorbar
% caxis([-200 200])

%%  CROP POINT CLOUD, PLOT NORMAL VECTORS, FIND SURFACE POINTS CORRESPONDING ONLY TO
% NORMALS THAT ARE PARALLEL TO CAMERA VIEWING VECTOR



%crop visible axial slice
points3d = pointCloud_dose_vis.Location;
points3d_1 = points3d(points3d(:, 3) > -200, :);
is = find(points3d(:,3)> -200);
pcCrop = pointCloud(points3d_1);
pcCrop.Intensity = pointCloud_dose_vis.Intensity(is);
v_iso_2D = world2img(pcCrop.Location,camExtrinsics,camIntrinsics);

ptCloud = pcCrop;
normals = pcnormals(ptCloud);
%normal origin points
x = ptCloud.Location(1:10:end,1);
y = ptCloud.Location(1:10:end,2);
z = ptCloud.Location(1:10:end,3);
% x = v_iso(1:10:end, 1);
% y = v_iso(1:10:end, 2);
% z = v_iso(1:10:end, 3);
% % directions
u = normals(1:10:end,1);
v = normals(1:10:end,2);
w = normals(1:10:end,3);
% u = N_ind(1:10:end, 1);
% v = N_ind(1:10:end, 2);
% w = N_ind(1:10:end, 3);

%we only want outward facing normals towards the camera
p1 = [-1615.4530911003208 -1068.4646556278083 -754.6003412378782]; %campos
p2 = [-1.7801874301633234 -17.336446178554752 27.650470285178017]; %camtarget, isocenter
cam_vect = p2-p1; %points FROM isocenter TO camera
norm_vs = [u, v,w];
cam_vects = -1.*cam_vect.*ones(length(norm_vs), 3);
sdot = sign(dot(norm_vs, cam_vects,3));
s_tot = sum(sdot,2);
pos_s = find(s_tot ==3);

x = x(pos_s);
y = y(pos_s);
z = z(pos_s);
u = u(pos_s);
v = v(pos_s);
w = w(pos_s);

% figure; pcshow(ptCloud,'BackgroundColor',[0 0 0],'ColorSource','Intensity'); hold on; 
% % Use external settings to set up camera
% camup([0.4675053886614372 -0.8623565422631594 0.19437053683222144]);
% campos([-1615.4530911003208 -1068.4646556278083 -754.6003412378782]);
% camtarget([-1.7801874301633234 -17.336446178554752 27.650470285178017]);
% camva(11.337405759636438)
% % hold on
% % quiver3(x,y,z,u,v,w, 20);
% 
% hold on; pcshow(p2, 'MarkerSize', 600)
% hold on; pcshow(p1, 'MarkerSize', 600)
% %hold on ; line([p2(1), p1(1)], [p2(2), p1(2)], [p2(3), p1(3)])
% vect =(p2 - p1)./norm(p2-p1);
% hold on; quiver3(p2(1),p2(2),p2(3),-1.*vect(1), -1.*vect(2), -1.*vect(3), 20)

% find only normals that are parallel to camera vector ----------------
% or find vectors within a certain
%cam vect direction 

inds = [];
crosses = ones(length(u), 3);
j = 0;
for i = 1:length(u) 
    v1 = [u(i), v(i), w(i)];
   % crosses(i, :) = cross(v1, cam_vect);
    a = atan2d(norm(cross(v1, -1.*cam_vect)),dot(v1, -1.*cam_vect)); % Angle in radians
  %  tf = isequal(cross(v1, cam_vect), [0 0 0]);
    if a <= 5 %if vectors are less than 5 degrees apart /or/ %if vectors are parallel 
        j = j+1;
        inds(j) = i;
    end
end
%

figure('Name','PointCloud','Position',[50 50 1200 1920]);
pcshow(pcCrop,'BackgroundColor',[0 0 0],'ColorSource','Intensity');
colormap(gray)
axis off 

% Use external settings to set up camera
camup([0.4675053886614372 -0.8623565422631594 0.19437053683222144]);
campos([-1615.4530911003208 -1068.4646556278083 -754.6003412378782]);
camtarget([-1.7801874301633234 -17.336446178554752 27.650470285178017]);
camva(11.337405759636438)

hold on 
x2 = x(inds);
y2 = y(inds);
z2 = z(inds);


quiver3(x(inds),y(inds),z(inds),u(inds),v(inds),w(inds), 10, 'LineWidth', 3);
hold off

    %%

%convert to image coords and plot visible point cloud

%v_iso_2D = world2img(ptCloud.Location,worldPose,camIntrinsics);

[C, ia, ib] = intersect( [x2, y2, z2], pointCloud_dose_vis.Location, 'rows' );
locs2 =pointCloud_dose_vis.Location;
ints2 = pointCloud_dose_vis.Intensity;
ints2(ib) = 2*max(pointCloud_dose_vis.Intensity);
ptCloud2 = pointCloud(locs2);
ptCloud2.Intensity = ints2;
%
v_iso_2D = world2img(pointCloud_dose_vis.Location,camExtrinsics,camIntrinsics);
%v_iso_2D = world2img(ptCloud2.Location,camExtrinsics,camIntrinsics);
figure;
%
scatter(v_iso_2D(:,1),v_iso_2D(:,2),4,ptCloud2.Intensity,'filled');

colormap(jet)
axis equal
set(gca,'Ydir','reverse')
xlim([1,1920]);
ylim([1,1200]);

% Interpolate visible point cloud to image matrix

F = scatteredInterpolant(v_iso_2D(:,1), v_iso_2D(:,2), pointCloud_dose_vis.Intensity,'linear','none');


rows = 1200;
columns = 1920;

[xGrid, yGrid] = meshgrid(1:columns, 1:rows);
xq = xGrid(:);
yq = yGrid(:);

vq = F(xq, yq);
ctImage_points = reshape(vq, rows, columns);
figure; imagesc(ctImage_points)
toc