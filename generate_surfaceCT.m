%% Generate Surface CT maps - Savannah Decker Fall 2023
% 
addpath('/Users/f004mg7/Library/CloudStorage/GoogleDrive-savannah.m.decker.th@dartmouth.edu/My Drive/Research/Scripts')
load('cameraCalib.mat') %camera parameters for Moffitt right cam

% direct to folder containing all patients of interest
cd('/Volumes/SavannahSSD/Moffitt/All Right Breast Supine')
%ctData2 = struct;

temp1 = pwd; files1 = dir(temp1);
PTfolders = {};
for f = 1:size(files1,1)
    if startsWith(files1(f).name, 'PT')
        PTfolders = [PTfolders;files1(f).name];
    end
end

for j =1:length(PTFolders)
    tic
    pt = (ctData4(j).Patient);
    cd((fullfile(temp1, pt, 'dicomData')))
    disp(strcat('Currently running: ', pt))
    disp(strcat('# of Patients Remaining: ', num2str(length(ctData4)-j)))

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

    planinfo = dicominfo(fullfile(temp,RPfile));
    fx = planinfo.FractionGroupSequence.Item_1.NumberOfFractionsPlanned;
    iso = planinfo.BeamSequence.Item_1.ControlPointSequence.Item_1.IsocenterPosition;


    % Read CT and header info -----------------00+

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
    try

 % if .tif file
    cd("ct/")
    ctvol_in = tiffreadVolume('ct.tif'); %some stacks need to be read in with FIJI and saved as a tiff stack. i.e.: test = tiffreadVolume('pt31ctstack.tif');
    ctimpospat = ctinfo.ImagePositionPatient; 
    ctres = [ctinfo.PixelSpacing', ctinfo.SliceThickness] ;
    ctvol = double(squeeze(ctvol_in));
% %
%  if using medicalVolume function
    % ctvol_in = medicalVolume(fullfile(temp, 'ct'));
    % ctimpospat = ctvol_in.VolumeGeometry.Position(1,:); %image position patient, x,y position of top left pixel
    % ctres= ctvol_in.VoxelSpacing;
    % ctvol = double(squeeze(dicomreadVolume(fullfile(temp, 'ct'))));

    ctvol = ctvol + offset;


    % % Read dose and header info -------------------
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

    ctvol = imresize3(ctvol, [size(ctvol,1), size(ctvol,2), 3*size(ctvol,3)]);
    ctres(3) = ctres(3)/3;
    dosesum = imresize3(dosesum, [size(dosesum,1), size(dosesum,2), 3*size(dosesum,3)]);
    doseres(3) = doseres(3)/3;
    xx = ctimpospat(1): ctres(1): ctimpospat(1)+ctres(1)*(size(ctvol,2)-1);
    yy = ctimpospat(2): ctres(2): ctimpospat(2)+ctres(2)*(size(ctvol,1)-1);
    zz = ctimpospat(3): ctres(3): ctimpospat(3)+ctres(3)*(size(ctvol,3)-1);
    xx = xx - iso(1);
    yy = yy - iso(2);
    zz = zz - iso(3) - 73.5; %try adding
    [ctXX,ctYY,ctZZ] = meshgrid(xx,yy,zz); 
    ctMesh.ctXX = ctXX;
    ctMesh.ctYY = ctYY;
    ctMesh.ctZZ = ctZZ;
    def_offset = -1000;
    ctvol_temp = ctvol;
    ctvol(400:end,:,:)= def_offset;

    ctvol(:, :, 1:end) = ctvol(:, :, end:-1:1); %needed if read in as .tif

    HUthres = -150;  
    [mm,nn,pp] = size(ctvol);
    [X,Y,Z] = meshgrid(1:nn, 1:mm, 1:pp); %index coords, interpolating to 1 mm
    surf_iso = isosurface(ctMesh.ctXX,ctMesh.ctYY,ctMesh.ctZZ,ctvol, HUthres); %iso coords
    surf_ind = isosurface(X,Y,Z,ctvol, HUthres); %index coords
    v_iso = surf_iso.vertices;
    v_ind = surf_ind.vertices;
    N_ind = isonormals(X,Y,Z,ctvol,v_ind, 'reverse');
 %  View isosurface
   % figure;
   % patch('Faces',surf_iso.faces,'Vertices',surf_iso.vertices,'FaceColor','flat','FaceVertexCData',[0 1 1],'EdgeColor','none');
   % set(gca,'DataAspectRatio',[1,1,1]);
   % light 
    ctColors = zeros(size(v_ind,1),1);
    f = waitbar(0,'Please wait...');
    Nn = numel(ctColors);
    d1 = 1; % Begin sampling at this location (in mm), assuming 1 mm spacing
    d2 =5; % End sampling at this location inside the patient
    for i=1:Nn
        coords = round(v_ind(i,:));
        xx = coords(1);
        yy = coords(2);
        zz = coords(3);
        vec = -1*N_ind(i,:);
        vec = vec/norm(vec); % Normalized inverse vector to follow for sampling
        sampCT = zeros(numel(d1:d2),1);
        k = 1;
        for d = d1:d2 % change sample depth, loop through sampling depth
            veci = round(d*vec);
            xxi = xx+veci(1);
            yyi = yy+veci(2);
            zzi = zz+veci(3);
            % To avoid going outside of volume
            if(xxi <= 0 || xxi > max(nn(:))) xxi = 1; end
            if(yyi <= 0 || yyi > max(mm(:)))  yyi = 1; end
            if(zzi <= 0 || zzi > max(pp(:))) zzi = 1; end
             sampCT(k) = ctvol(yyi,xxi,zzi);
             k = k+1;
        end
        ctColors(i)     =   mean(sampCT(:),'omitnan');
        if mod(i,5000) == 0
            perc = i/Nn;
            waitbar(perc,f,['Please wait (',num2str(100*perc,2),'%) ...'])
        end
    end
    close(f)
    % Color the faces with the surface dose or HU values-------------------------
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

    %Moffitt Right Cam
    camup([0.4675053886614372 -0.8623565422631594 0.19437053683222144]);
    campos([-1615.4530911003208 -1068.4646556278083 -754.6003412378782]);
    camtarget([-1.7801874301633234 -17.336446178554752 27.650470285178017]);
    camva(11.337405759636438)
    camproj('perspective')

    camup('manual')
    clim([-200 100])
    colormap(gray)
    ax = gca; fig = gcf; ax.Position = [0 0 1 1]; fig.Position = [0 0 1528 955];
    frame = getframe(gcf);
    imwrite(frame.cdata,'test.png')
    close

    test_im = sum(imread('test.png'),3);
    test_im_rescaled = (rescale(test_im, -200, 100));
    figure
    imagesc(test_im_rescaled); axis image; 
    clim([-200 100])
    colormap(gray)
    % 
    % pointCloud_dose = pointCloud(v_iso);
    % pointCloud_dose.Intensity = ctColors;
    % [pointCloud_dose_vis,visInd] = removeHiddenPoints(pointCloud_dose,...
    %     [-1615.4530911003208 -1068.4646556278083 -754.6003412378782], RadiusScale=5);
    % figure('Name','PointCloud','Position',[50 50 1200 1920]);
    % pcshow(pointCloud_dose_vis, 'BackgroundColor',[0 0 0],'ColorSource','Intensity');
    % colormap(gray)
    % axis off 
    % % % Use external settings to set up camera
    % camup([0.4675053886614372 -0.8623565422631594 0.19437053683222144]);
    % campos([-1615.4530911003208 -1068.4646556278083 -754.6003412378782]);
    % camtarget([-1.7801874301633234 -17.336446178554752 27.650470285178017]);
    % camva(11.337405759636438)
    % v_iso_2D = world2img(pointCloud_dose_vis.Location,camExtrinsics,camIntrinsics);
    % F = scatteredInterpolant(v_iso_2D(:,1), v_iso_2D(:,2), pointCloud_dose_vis.Intensity,'linear','none');
    % rows = 1200;
    % columns = 1920;
    % [xGrid, yGrid] = meshgrid(1:columns, 1:rows);
    % xq = xGrid(:);
    % yq = yGrid(:);ces
    % vq = F(xq, yq);
    % ctImage = reshape(vq, rows, columns);
    % ctData2(j).Patient = pt;
    % ctData2(j).ctvol = ctvol;
     ctData4(j).ctImage =test_im_rescaled;
     close
    toc
   catch
       continue
   end

end
% cd('/Users/f004mg7/Downloads/Moffitt/All Right Breast Supine')
% save('ctData2.mat', 'ctData2')

%% view all
% figure; tiledlayout(3,4, 'TileSpacing', 'Compact')
% for i = [7,11,16,17,26,27,29,30,31,33,34,36]
%     nexttile
%     imagesc(PTfolders{i,2})
%     colormap gray; caxis([-200 100]); axis off
%     title(PTfolders{i,1})
% end
% %% view new data
% figure; tiledlayout(3,5, 'TileSpacing', 'Compact')
% for i = [8, 9, 11, 13, 14]%[1, 3, 4,7, 12, 15, 2, 5,10,17,19, 21,25]
%     nexttile 
%     imagesc(ctData3(i).ctImage)
%     colormap gray; clim([-200 100]); axis off
%     title(ctData3(i).Patient)
% end
%% remove tape data

for i = 1:length(ctData4)
%    ctData3(i).ctImage = ctData1(i).ctImage ;
    temp=ctData4(i).ctImage;
    temp(temp==100)=NaN;
     figure; imagesc(temp); colormap gray; clim([-200 100]); axis off
     title(ctData4(i).Patient)
     ctData4(i).CTnotape = temp;
  % %  close
end
%% add new data to exisiting struct
js = [7,11,16,17,24,26,27,29,30,31,33:37];
for i = 1:length(js)
    ind = 27+i;
    ctData4(ind).Patient = cell2mat(PTfolders(js(i), 1));
    ctData4(ind).ctImage = cell2mat(PTfolders(js(i), 2));
    ctData4(ind).CTnotape = cell2mat(PTfolders(js(i), 3));
end
