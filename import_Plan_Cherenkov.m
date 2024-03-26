%import Cherenkov beam acqs and indivudal beam plan images from CDose

addpath('/Users/f004mg7/Library/CloudStorage/GoogleDrive-savannah.m.decker.th@dartmouth.edu/My Drive/Research/Scripts')
load('cameraCalib.mat') %camera parameters for Moffitt right cam

% direct to folder containing all patients of interest
cd('/Volumes/SavannahSSD/Moffitt/All Right Breast Supine')
load('plan2chkv_reg.mat')

temp1 = pwd; files1 = dir(temp1);
PTfolders = {};
for f = 1:size(files1,1)
    if startsWith(files1(f).name, 'PT')
        PTfolders = [PTfolders;files1(f).name];
    end
end

%% Cherenkov

for j = 1:numel(PTfolders)
    pt = cell2mat(PTfolders(j));
    %import Cherenkov
    cd((fullfile(temp1, pt, 'Cherenkov')))
    % disp(strcat('Currently running: ', pt))
    % disp(strcat('# of Patients Remaining: ', num2str(length(PTfolders)-j)))
    acqs = dir(pwd); 
    for k = 1:length(acqs)
        if startsWith(acqs(k).name, '20') ==1
            cd(acqs(k).name)
            acq_num = acqs(k).name(end-3:end);
            ctData5(j).(acq_num) = read_dovi_c('meas_cam1.dovi', 's1');
            cd ..
        end
    end
    cd ..
end
%% Dose
cd('/Volumes/SavannahSSD/Moffitt/All Right Breast Supine')
for j = [12]
    pt = (ctData5(j).Patient);
    %import Cherenkov
    cd((fullfile(temp1, pt, 'planDose')))
    beams = dir(pwd); 
    for k = 1:length(beams)
        if startsWith(beams(k).name, 'be') ==1
            cd(beams(k).name)
            doseMap = zeros(1200 ,1920);
            temp = double(rgb2gray(imread('screenshot.png')));
            figure; imagesc(temp);
            h=imrect(); m = createMask(h);
            caxisMax = input('Max Colorbar Value:');
            temp2 = rescale(temp,0, caxisMax) .*imcomplement(m);
            temp3 = imwarp(temp2, imref2d(size(plan_check)), movingReg.Transformation, 'OutputView', movingReg.SpatialRefObj, 'SmoothEdges', true);
            doseMap = doseMap + temp3;
            beam = beams(k).name;
            ctData5(j).(beam) = doseMap;
            cd ..
        end
        clear doseMap temp temp2 temp
    end
    cd ..
end
      

%% show Cherenkov and dose maps
%% PT 1 %%%
figure()
tiledlayout(3,4, 'TileSpacing','tight')
j = 1;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ1); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam D')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ2); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam D1')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ6); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam E')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ7);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam E1')); set(gca,'FontSize', 18)
%Dose
nexttile
imagesc(ctData5(j).beamD); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam D - 6x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).beamD1); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam D1 - 15x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).beamE);  colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam E - 15x RPO')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).beamE1); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam E1 - 6x RPO')); set(gca,'FontSize', 18)

% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ1); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (ctData5(j).beamD .* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam D Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ2); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (ctData5(j).beamD1 .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam D1 Norm')); set(gca,'FontSize', 18)
nexttile
bn3 = double(ctData5(j).ACQ6); msk=double(imbinarize(bn3, 0.2.*prctile(bn3(:), 99.9)));
bn3 = (bn3.*msk) ./ (ctData5(j).beamE .* msk);
imagesc(bn3);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam E Norm')); set(gca,'FontSize', 18)
nexttile
bn4 = double(ctData5(j).ACQ7); msk=double(imbinarize(bn4, 0.2.*prctile(bn4(:), 99.9)));
bn4 = (bn4.*msk) ./ (ctData5(j).beamE1 .* msk);
imagesc(bn4);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam E1 Norm')); set(gca,'FontSize', 18)

% save dose norms
ctData5(j).RPO6Xdn = bn4;
ctData5(j).RPO15Xdn = bn3;
ctData5(j).LAO6Xdn = bn1;
ctData5(j).LAO15Xdn = bn2;

j=1;
ctData5(j).total6chkv = ctData5(j).ACQ1 + ctData5(j).ACQ7;
ctData5(j).total6dose = ctData5(j).beamD + ctData5(j).beamE1;

ctData5(j).total15chkv = ctData5(j).ACQ2 + ctData5(j).ACQ6;
ctData5(j).total15dose = ctData5(j).beamD1 + ctData5(j).beamE;
%% PT 3 %%% error found in spreadsheet...talk to Diego
figure
tiledlayout(3,2, 'TileSpacing','tight')
j = 7;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ1); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam E')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ5);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam D')); set(gca,'FontSize', 18)
%Dose
nexttile
imagesc(medfilt2(ctData5(j).beamE, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam E - 15x RPO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamD, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam D - 6x LAO')); set(gca,'FontSize', 18)

% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ1); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (medfilt2(ctData5(j).beamE, [30, 30]).* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam E Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ5); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (medfilt2(ctData5(j).beamD, [30, 30]) .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam D Norm')); set(gca,'FontSize', 18)

ctData5(j).RPO6Xdn = [];
ctData5(j).RPO15Xdn = bn1;
ctData5(j).LAO6Xdn = bn2;
ctData5(j).LAO15Xdn = [];
%% PT 8 %%%
figure()
tiledlayout(3,3, 'TileSpacing','tight')
j = 29;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ1); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam A')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ2); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam B')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ5); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam C')); set(gca,'FontSize', 18)

%Dose
nexttile
imagesc(ctData5(j).beamA); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam A - 6x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).beamB); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam B - 15x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).beamC);  colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam C - 15x RPO')); set(gca,'FontSize', 18)


% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ1); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (ctData5(j).beamA .* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam A Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ2); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (ctData5(j).beamB .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam B Norm')); set(gca,'FontSize', 18)
nexttile
bn3 = double(ctData5(j).ACQ5); msk=double(imbinarize(bn3, 0.2.*prctile(bn3(:), 99.9)));
bn3 = (bn3.*msk) ./ (ctData5(j).beamC .* msk);
imagesc(bn3);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam C Norm')); set(gca,'FontSize', 18)


% save dose norms
ctData5(j).RPO6Xdn = [];
ctData5(j).RPO15Xdn = bn3;
ctData5(j).LAO6Xdn = bn1;
ctData5(j).LAO15Xdn = bn2;


ctData5(j).total6chkv = ctData5(j).ACQ1;
ctData5(j).total6dose = ctData5(j).beamA;

ctData5(j).total15chkv = ctData5(j).ACQ2 + ctData5(j).ACQ5;
ctData5(j).total15dose = ctData5(j).beamB + ctData5(j).beamC;
%% PT 12 %%% 
figure
tiledlayout(3,2, 'TileSpacing','tight')
j = 2;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ1); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam C1')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ2);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam D1')); set(gca,'FontSize', 18)
%Dose
nexttile
imagesc(medfilt2(ctData5(j).beamC1, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam C1 - 15x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamD1, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam D1 - 15x RPO')); set(gca,'FontSize', 18)

% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ3); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (medfilt2(ctData5(j).beamA, [30, 30]).* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam A Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ4); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (medfilt2(ctData5(j).beamB, [30, 30]) .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam B Norm')); set(gca,'FontSize', 18)

ctData5(j).RPO6Xdn = [];
ctData5(j).RPO15Xdn = bn2;
ctData5(j).LAO6Xdn = [];
ctData5(j).LAO15Xdn = bn1;
%% PT 14 %%% 
figure
tiledlayout(3,3, 'TileSpacing','tight')
j = 2;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ3); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam A')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ4);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam A1')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ5);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam B')); set(gca,'FontSize', 18)
%Dose
nexttile
imagesc(medfilt2(ctData5(j).beamA, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam A - 6x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamA1, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam A1 - 15x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamB, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam B - 15x RPO')); set(gca,'FontSize', 18)

% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ3); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (medfilt2(ctData5(j).beamA, [30, 30]).* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam A Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ4); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (medfilt2(ctData5(j).beamA1, [30, 30]) .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam A1 Norm')); set(gca,'FontSize', 18)
nexttile
bn3 = double(ctData5(j).ACQ5); msk=double(imbinarize(bn3, 0.2.*prctile(bn3(:), 99.9)));
bn3 = (bn3.*msk) ./ (medfilt2(ctData5(j).beamB, [30, 30]) .* msk);
imagesc(bn3);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam B Norm')); set(gca,'FontSize', 18)

ctData5(j).RPO6Xdn = [];
ctData5(j).RPO15Xdn = bn3;
ctData5(j).LAO6Xdn = bn1;
ctData5(j).LAO15Xdn = bn2;
%
j=3;
ctData5(j).total6chkv = ctData5(j).ACQ3 ;
ctData5(j).total6dose = ctData5(j).beamA;

ctData5(j).total15chkv = ctData5(j).ACQ4 + ctData5(j).ACQ5;
ctData5(j).total15dose = ctData5(j).beamA1 + ctData5(j).beamB;

%% PT 15 %%% 
figure
tiledlayout(3,2, 'TileSpacing','tight')
j = 4;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ3); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam A')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ4));  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam B')); set(gca,'FontSize', 18)
%Dose
nexttile
imagesc(medfilt2(ctData5(j).beamA, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam A - 15x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamB, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam B - 15x RPO')); set(gca,'FontSize', 18)

% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ3); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (medfilt2(ctData5(j).beamA, [30, 30]).* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam A Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ4); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (medfilt2(ctData5(j).beamB, [30, 30]) .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam B Norm')); set(gca,'FontSize', 18)

ctData5(j).RPO6Xdn = [];
ctData5(j).RPO15Xdn = bn2;
ctData5(j).LAO6Xdn = [];
ctData5(j).LAO15Xdn = bn1;
%
j=4;
ctData5(j).total15chkv = ctData5(j).ACQ3 + ctData5(j).ACQ4;
ctData5(j).total15dose = ctData5(j).beamA + ctData5(j).beamB;

%% PT 17 %%% 
figure
tiledlayout(3,2, 'TileSpacing','tight')
j = 5;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ2); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam D')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ6);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam C')); set(gca,'FontSize', 18)
%Dose
nexttile
imagesc(medfilt2(ctData5(j).beamD, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam D - 15x RPO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamC, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam C - 6x LAO')); set(gca,'FontSize', 18)

% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ2); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (medfilt2(ctData5(j).beamD, [30, 30]).* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam D Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ6); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (medfilt2(ctData5(j).beamC, [30, 30]) .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam C Norm')); set(gca,'FontSize', 18)

ctData5(j).RPO6Xdn = [];
ctData5(j).RPO15Xdn = bn1;
ctData5(j).LAO6Xdn = bn2;
ctData5(j).LAO15Xdn = [];

j=5;
ctData5(j).total6chkv = ctData5(j).ACQ6;
ctData5(j).total6dose = ctData5(j).beamC;

ctData5(j).total15chkv = ctData5(j).ACQ2;
ctData5(j).total15dose = ctData5(j).beamD;
%% PT 25 %%% 
figure
tiledlayout(3,3, 'TileSpacing','tight')
j = 6;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ1); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam C')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ3);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam C1')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ6);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam D')); set(gca,'FontSize', 18)
%Dose
nexttile
imagesc(medfilt2(ctData5(j).beamC, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam C - 6x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamC1, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam C1 - 15x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamD, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam C1 - 15x RPO')); set(gca,'FontSize', 18)

% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ1); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (medfilt2(ctData5(j).beamC, [30, 30]).* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam C Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ3); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (medfilt2(ctData5(j).beamC1, [30, 30]) .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam C1 Norm')); set(gca,'FontSize', 18)
nexttile
bn3 = double(ctData5(j).ACQ6); msk=double(imbinarize(bn3, 0.2.*prctile(bn3(:), 99.9)));
bn3 = (bn3.*msk) ./ (medfilt2(ctData5(j).beamD, [30, 30]) .* msk);
imagesc(bn3);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam D Norm')); set(gca,'FontSize', 18)

ctData5(j).RPO6Xdn = [];
ctData5(j).RPO15Xdn = bn3;
ctData5(j).LAO6Xdn = bn1;
ctData5(j).LAO15Xdn= bn2;

j = 6;
ctData5(j).total15chkv = ctData5(j).ACQ3 + ctData5(j).ACQ6;
ctData5(j).total15dose = ctData5(j).beamC1 + ctData5(j).beamD;
ctData5(j).total6chkv = ctData5(j).ACQ1;
ctData5(j).total6dose = ctData5(j).beamC;
%% PT 30 %%% 
figure
tiledlayout(3,2, 'TileSpacing','tight')
j = 8;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ2); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam A')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ5);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam B')); set(gca,'FontSize', 18)
%Dose
nexttile
imagesc(medfilt2(ctData5(j).beamA, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam A - 15x RPO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamB, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam B - 15x LAO')); set(gca,'FontSize', 18)

% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ2); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (medfilt2(ctData5(j).beamA, [30, 30]).* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam A Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ5); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (medfilt2(ctData5(j).beamB, [30, 30]) .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam B Norm')); set(gca,'FontSize', 18)

ctData5(j).RPO6Xdn = [];
ctData5(j).RPO15Xdn = bn1;
ctData5(j).LAO6Xdn = [];
ctData5(j).LAO15Xdn = bn2;

j = 8;
ctData5(j).total15chkv = ctData5(j).ACQ2 + ctData5(j).ACQ5;
ctData5(j).total15dose = ctData5(j).beamA + ctData5(j).beamB;
%% PT 26 %%%
figure()
tiledlayout(3,3, 'TileSpacing','tight')
j = 4;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ1); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam A')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ2); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam B')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ3); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam C')); set(gca,'FontSize', 18)

%Dose
nexttile
imagesc(ctData5(j).beamA); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam A - 6x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).beamB); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam B - 6x RPO')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).beamC);  colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam C - 15x RPO')); set(gca,'FontSize', 18)


% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ1); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (ctData5(j).beamA .* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam A Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ2); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (ctData5(j).beamB .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam B Norm')); set(gca,'FontSize', 18)
nexttile
bn3 = double(ctData5(j).ACQ3); msk=double(imbinarize(bn3, 0.2.*prctile(bn3(:), 99.9)));
bn3 = (bn3.*msk) ./ (ctData5(j).beamC .* msk);
imagesc(bn3);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam C Norm')); set(gca,'FontSize', 18)


% save dose norms
ctData5(j).RPO6Xdn = bn2;
ctData5(j).RPO15Xdn = bn3;
ctData5(j).LAO6Xdn = bn1;
ctData5(j).LAO15Xdn = [];


ctData5(j).total6chkv = ctData5(j).ACQ1 +ctData5(j).ACQ2;
ctData5(j).total6dose = ctData5(j).beamA +ctData5(j).beamB;

ctData5(j).total15chkv = ctData5(j).ACQ3;
ctData5(j).total15dose = ctData5(j).beamC;
% PT 31 %%% 
figure
tiledlayout(3,2, 'TileSpacing','tight')
j = 9;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ2); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam A')); set(gca,'FontSize', 18)
%Dose
nexttile
imagesc(medfilt2(ctData5(j).beamA, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam A - 15x RPO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamB, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam B - 15x LAO')); set(gca,'FontSize', 18)

% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ2); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (medfilt2(ctData5(j).beamA, [30, 30]).* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam A Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ5); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (medfilt2(ctData5(j).beamB, [30, 30]) .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam B Norm')); set(gca,'FontSize', 18)

ctData5(j).RPO6Xdn = [];
ctData5(j).RPO15Xdn = bn1;
ctData5(j).LAO6Xdn = [];
ctData5(j).LAO15Xdn = bn2;

j = 9;
ctData5(j).total6chkv = ctData5(j).ACQ1+ ctData5(j).ACQ2+ ctData5(j).ACQ4;
ctData5(j).total6dose = ctData5(j).beamA;
%% PT 34 %%%
figure()
tiledlayout(3,1, 'TileSpacing','tight')
j = 6;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ4); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam A')); set(gca,'FontSize', 18)

%Dose
nexttile
imagesc(ctData5(j).beamA); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam A - 6x RPO')); set(gca,'FontSize', 18)

% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ4); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (ctData5(j).beamA .* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam A Norm')); set(gca,'FontSize', 18)


% save dose norms
ctData5(j).RPO6Xdn = bn1;



ctData5(j).total6chkv = ctData5(j).ACQ4;
ctData5(j).total6dose = ctData5(j).beamA ;


%% PT 35 %%%
figure
tiledlayout(3,2, 'TileSpacing','tight')
j = 7;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ1); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam A')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ2);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam B')); set(gca,'FontSize', 18)
%Dose
nexttile
imagesc(medfilt2(ctData5(j).beamA, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam A - 6x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamB, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam B - 6x RPO')); set(gca,'FontSize', 18)

% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ1); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (medfilt2(ctData5(j).beamA, [30, 30]).* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam A Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ2); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (medfilt2(ctData5(j).beamB, [30, 30]) .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam B Norm')); set(gca,'FontSize', 18)

ctData5(j).RPO6Xdn = bn2;
ctData5(j).RPO15Xdn = [];
ctData5(j).LAO6Xdn = bn1;
ctData5(j).LAO15Xdn = [];

j = 10;
ctData5(j).total6chkv = ctData5(j).ACQ1+ ctData5(j).ACQ2;
ctData5(j).total6dose = ctData5(j).beamA + ctData5(j).beamB;

%% PT 38 %%%
figure
tiledlayout(3,2, 'TileSpacing','tight')
j = 11;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ5); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam F')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ6);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam E')); set(gca,'FontSize', 18)
%Dose
nexttile
imagesc(medfilt2(ctData5(j).beamF, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam A - 6x IMRT')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamE, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam B - 6x IMRT')); set(gca,'FontSize', 18)

% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ5); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (medfilt2(ctData5(j).beamF, [30, 30]).* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam F Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ6); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (medfilt2(ctData5(j).beamE, [30, 30]) .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam E Norm')); set(gca,'FontSize', 18)

% ctData5(j).RPO6Xdn = bn2;
% ctData5(j).RPO15Xdn = [];
% ctData5(j).LAO6Xdn = bn1;
% ctData5(j).LAO15Xdn = [];

j = 11;
ctData5(j).total6chkv = ctData5(j).ACQ5+ ctData5(j).ACQ6;
ctData5(j).total6dose = ctData5(j).beamE + ctData5(j).beamF;
%% PT 39 %%% 
figure
tiledlayout(3,3, 'TileSpacing','tight')
j = 12;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ5); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam C')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ6);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam C1')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ7);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam D')); set(gca,'FontSize', 18)
%Dose
nexttile
imagesc(medfilt2(ctData5(j).beamC, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam C - 6x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamC1, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam C1 - 15x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamD, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam D - 15x RPO')); set(gca,'FontSize', 18)

% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ5); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (medfilt2(ctData5(j).beamC, [30, 30]).* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam C Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ6); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (medfilt2(ctData5(j).beamC1, [30, 30]) .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam C1 Norm')); set(gca,'FontSize', 18)
nexttile
bn3 = double(ctData5(j).ACQ7); msk=double(imbinarize(bn3, 0.2.*prctile(bn3(:), 99.9)));
bn3 = (bn3.*msk) ./ (medfilt2(ctData5(j).beamD, [30, 30]) .* msk);
imagesc(bn3);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam D Norm')); set(gca,'FontSize', 18)

ctData5(j).RPO6Xdn = [];
ctData5(j).RPO15Xdn = bn3;
ctData5(j).LAO6Xdn = bn1;
ctData5(j).LAO15Xdn = bn2;

j = 12;
ctData5(j).total15chkv = ctData5(j).ACQ6 + ctData5(j).ACQ7;
ctData5(j).total15dose = ctData5(j).beamC1 + ctData5(j).beamD;
ctData5(j).total6chkv = ctData5(j).ACQ5;
ctData5(j).total6dose = ctData5(j).beamC;
%% PT 41 %%%
figure
tiledlayout(3,2, 'TileSpacing','tight')
j = 13;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ2); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam A')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ4);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam C')); set(gca,'FontSize', 18)
%Dose
nexttile
imagesc(medfilt2(ctData5(j).beamA, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam A - 15X RPO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamC, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam C - 15X LAO')); set(gca,'FontSize', 18)

% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ2); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (medfilt2(ctData5(j).beamA, [30, 30]).* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam A Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ4); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (medfilt2(ctData5(j).beamC, [30, 30]) .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam C Norm')); set(gca,'FontSize', 18)

ctData5(j).RPO6Xdn =[];
ctData5(j).RPO15Xdn = bn1;
ctData5(j).LAO6Xdn = [];
ctData5(j).LAO15Xdn = bn2;

j = 13;
ctData5(j).total15chkv = ctData5(j).ACQ2+ ctData5(j).ACQ4;
ctData5(j).total15dose = ctData5(j).beamA + ctData5(j).beamC;
%% PT 47 %%%
figure()
tiledlayout(3,3, 'TileSpacing','tight')
j = 13;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ1); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam A')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ3); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam B')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ2); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam C')); set(gca,'FontSize', 18)

%Dose
nexttile
imagesc(ctData5(j).beamA); colormap(cmap1); axis off; clim([0 350])
title(strcat(pt, ' Beam A - 6x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).beamB); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam B - 6x RPO')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).beamC);  colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam C - 15x RPO')); set(gca,'FontSize', 18)


% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ1); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (ctData5(j).beamA .* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam A Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ3); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (ctData5(j).beamB .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam B Norm')); set(gca,'FontSize', 18)
nexttile
bn3 = double(ctData5(j).ACQ2); msk=double(imbinarize(bn3, 0.2.*prctile(bn3(:), 99.9)));
bn3 = (bn3.*msk) ./ (ctData5(j).beamC .* msk);
imagesc(bn3);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam C Norm')); set(gca,'FontSize', 18)


% save dose norms
ctData5(j).RPO6Xdn = bn2;
ctData5(j).RPO15Xdn = bn3;
ctData5(j).LAO6Xdn = bn1;
ctData5(j).LAO15Xdn = [];


ctData5(j).total6chkv = ctData5(j).ACQ3 +ctData5(j).ACQ1;
ctData5(j).total6dose = ctData5(j).beamA +ctData5(j).beamB;

ctData5(j).total15chkv = ctData5(j).ACQ2;
ctData5(j).total15dose = ctData5(j).beamC;
%% PT 49 %%% 
figure
tiledlayout(3,3, 'TileSpacing','tight')
j = 14;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ1); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam A')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ2);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam B')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ4);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam D')); set(gca,'FontSize', 18)
%Dose
nexttile
imagesc(medfilt2(ctData5(j).beamA, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam A - 6x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamB, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam B - 15x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamD, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam D - 15x RPO')); set(gca,'FontSize', 18)

% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ1); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (medfilt2(ctData5(j).beamA, [30, 30]).* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam A Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ2); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (medfilt2(ctData5(j).beamB, [30, 30]) .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam B Norm')); set(gca,'FontSize', 18)
nexttile
bn3 = double(ctData5(j).ACQ4); msk=double(imbinarize(bn3, 0.2.*prctile(bn3(:), 99.9)));
bn3 = (bn3.*msk) ./ (medfilt2(ctData5(j).beamD, [30, 30]) .* msk);
imagesc(bn3);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam D Norm')); set(gca,'FontSize', 18)

ctData5(j).RPO6Xdn = [];
ctData5(j).RPO15Xdn = bn3;
ctData5(j).LAO6Xdn = bn1;
ctData5(j).LAO15Xdn= bn2;

j = 14;
ctData5(j).total15chkv = ctData5(j).ACQ2 + ctData5(j).ACQ4;
ctData5(j).total15dose = ctData5(j).beamB + ctData5(j).beamD;
ctData5(j).total6chkv = ctData5(j).ACQ1;
ctData5(j).total6dose = ctData5(j).beamA;
%% PT 51 %%% 
figure
tiledlayout(3,3, 'TileSpacing','tight')
j = 15;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ1); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam C')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ4);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam D')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ5);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam D1')); set(gca,'FontSize', 18)
%Dose
nexttile
imagesc(medfilt2(ctData5(j).beamC, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam C - 6x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamD, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam D - 15x RPO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamD1, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam D1 - 6x RPO')); set(gca,'FontSize', 18)

% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ1); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (medfilt2(ctData5(j).beamC, [30, 30]).* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam C Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ4); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (medfilt2(ctData5(j).beamD, [30, 30]) .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam D Norm')); set(gca,'FontSize', 18)
nexttile
bn3 = double(ctData5(j).ACQ5); msk=double(imbinarize(bn3, 0.2.*prctile(bn3(:), 99.9)));
bn3 = (bn3.*msk) ./ (medfilt2(ctData5(j).beamD1, [30, 30]) .* msk);
imagesc(bn3);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam D1 Norm')); set(gca,'FontSize', 18)

ctData5(j).RPO6Xdn = bn3;
ctData5(j).RPO15Xdn = bn2;
ctData5(j).LAO6Xdn = bn1;
ctData5(j).LAO15Xdn = [];

%% PT 53 %%% 
figure
tiledlayout(3,4, 'TileSpacing','tight')
j = 11;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ3); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam A')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ4);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam A1')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ1);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam B')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ2);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam B1')); set(gca,'FontSize', 18)
%Dose
nexttile
imagesc(medfilt2(ctData5(j).beamA, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam A - 6x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamA1, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam A1 - 15x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamB, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam B - 6x RPO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamB1, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam B1 - 15x RPO')); set(gca,'FontSize', 18)

% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ3); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (medfilt2(ctData5(j).beamA, [30, 30]).* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam A Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ4); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (medfilt2(ctData5(j).beamA1, [30, 30]) .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam A1 Norm')); set(gca,'FontSize', 18)
nexttile
bn3 = double(ctData5(j).ACQ1); msk=double(imbinarize(bn3, 0.2.*prctile(bn3(:), 99.9)));
bn3 = (bn3.*msk) ./ (medfilt2(ctData5(j).beamB, [30, 30]) .* msk);
imagesc(bn3);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam B Norm')); set(gca,'FontSize', 18)
nexttile
bn4 = double(ctData5(j).ACQ2); msk=double(imbinarize(bn4, 0.2.*prctile(bn4(:), 99.9)));
bn4 = (bn4.*msk) ./ (medfilt2(ctData5(j).beamB1, [30, 30]) .* msk);
imagesc(bn4);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam B1 Norm')); set(gca,'FontSize', 18)

ctData5(j).RPO6Xdn = bn3;
ctData5(j).RPO15Xdn = bn4;
ctData5(j).LAO6Xdn = bn1;
ctData5(j).LAO15Xdn = bn2;

j=17;
ctData5(j).total6chkv = ctData5(j).ACQ3 + ctData5(j).ACQ1;
ctData5(j).total6dose = ctData5(j).beamA + ctData5(j).beamB;

ctData5(j).total15chkv = ctData5(j).ACQ4 + ctData5(j).ACQ2;
ctData5(j).total15dose = ctData5(j).beamA1 + ctData5(j).beamB1;
%% PT 54 %%% 
figure
tiledlayout(3,3, 'TileSpacing','tight')
j = 12;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ3); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam A')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ4);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam B')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ7);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam D')); set(gca,'FontSize', 18)

%Dose
nexttile
imagesc(medfilt2(ctData5(j).beamA, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam A - 6x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamB, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam B - 15x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamD, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam D - 15x RPO')); set(gca,'FontSize', 18)

% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ3); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (medfilt2(ctData5(j).beamA, [30, 30]).* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam A Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ4); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (medfilt2(ctData5(j).beamB, [30, 30]) .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam B Norm')); set(gca,'FontSize', 18)
nexttile
bn3 = double(ctData5(j).ACQ7); msk=double(imbinarize(bn3, 0.2.*prctile(bn3(:), 99.9)));
bn3 = (bn3.*msk) ./ (medfilt2(ctData5(j).beamD, [30, 30]) .* msk);
imagesc(bn3);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam D Norm')); set(gca,'FontSize', 18)


ctData5(j).RPO6Xdn = [];
ctData5(j).RPO15Xdn = bn3;
ctData5(j).LAO6Xdn = bn1;
ctData5(j).LAO15Xdn = bn2;

j=12;
ctData5(j).total6chkv = ctData5(j).ACQ3 ;
ctData5(j).total6dose = ctData5(j).beamA ;

ctData5(j).total15chkv = ctData5(j).ACQ4 + ctData5(j).ACQ7;
ctData5(j).total15dose = ctData5(j).beamB + ctData5(j).beamD;

%% PT 55 %%% 
figure
tiledlayout(3,4, 'TileSpacing','tight')
j = 13;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ5); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam A')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ6);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam B')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ2);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam C')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ1);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam D')); set(gca,'FontSize', 18)
%Dose
nexttile
imagesc(medfilt2(ctData5(j).beamA, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam A - 6x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamB, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam B - 15x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamC, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam C - 6x RPO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamD, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam D - 15x RPO')); set(gca,'FontSize', 18)

% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ5); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (medfilt2(ctData5(j).beamA, [30, 30]).* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam A Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ6); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (medfilt2(ctData5(j).beamB, [30, 30]) .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam A1 Norm')); set(gca,'FontSize', 18)
nexttile
bn3 = double(ctData5(j).ACQ2); msk=double(imbinarize(bn3, 0.2.*prctile(bn3(:), 99.9)));
bn3 = (bn3.*msk) ./ (medfilt2(ctData5(j).beamC, [30, 30]) .* msk);
imagesc(bn3);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam B Norm')); set(gca,'FontSize', 18)
nexttile
bn4 = double(ctData5(j).ACQ1); msk=double(imbinarize(bn4, 0.2.*prctile(bn4(:), 99.9)));
bn4 = (bn4.*msk) ./ (medfilt2(ctData5(j).beamD, [30, 30]) .* msk);
imagesc(bn4);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam B Norm')); set(gca,'FontSize', 18)

ctData5(j).RPO6Xdn = bn3;
ctData5(j).RPO15Xdn = bn4;
ctData5(j).LAO6Xdn = bn1;
ctData5(j).LAO15Xdn = bn2;

j=13;
ctData5(j).total6chkv = ctData5(j).ACQ5 + ctData5(j).ACQ2;
ctData5(j).total6dose = ctData5(j).beamA + ctData5(j).beamC;

ctData5(j).total15chkv = ctData5(j).ACQ6 + ctData5(j).ACQ1;
ctData5(j).total15dose = ctData5(j).beamB + ctData5(j).beamD;
%% PT 56 %%% 
figure
tiledlayout(3,5, 'TileSpacing','tight')
j = 20;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ7); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam A')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ6);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam B')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ3);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam E')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ2);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam F')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ1);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam G')); set(gca,'FontSize', 18)
%Dose
nexttile
imagesc(medfilt2(ctData5(j).beamA, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam A - 15x RPO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamB, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam B - 15x RPO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamE, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam E - 6x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamF, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam F - 6x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamG, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam G - 15x LAO')); set(gca,'FontSize', 18)


% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ7); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (medfilt2(ctData5(j).beamA, [30, 30]).* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam A Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ6); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (medfilt2(ctData5(j).beamB, [30, 30]) .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam B Norm')); set(gca,'FontSize', 18)
nexttile
bn3 = double(ctData5(j).ACQ3); msk=double(imbinarize(bn3, 0.2.*prctile(bn3(:), 99.9)));
bn3 = (bn3.*msk) ./ (medfilt2(ctData5(j).beamE, [30, 30]) .* msk);
imagesc(bn3);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam E Norm')); set(gca,'FontSize', 18)
nexttile
bn4 = double(ctData5(j).ACQ2); msk=double(imbinarize(bn4, 0.2.*prctile(bn4(:), 99.9)));
bn4 = (bn4.*msk) ./ (medfilt2(ctData5(j).beamF, [30, 30]) .* msk);
imagesc(bn4);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam F Norm')); set(gca,'FontSize', 18)
nexttile
bn5 = double(ctData5(j).ACQ1); msk=double(imbinarize(bn5, 0.2.*prctile(bn5(:), 99.9)));
bn5 = (bn5.*msk) ./ (medfilt2(ctData5(j).beamG, [30, 30]) .* msk);
imagesc(bn5);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam G Norm')); set(gca,'FontSize', 18)

ctData5(j).RPO6Xdn = [];
ctData5(j).RPO15Xdn = bn1 + bn2;
ctData5(j).LAO6Xdn = bn3 + bn4;
ctData5(j).LAO15Xdn = bn5;


ctData5(j).total6chkv = ctData5(j).ACQ3 + ctData5(j).ACQ2;
ctData5(j).total6dose = ctData5(j).beamE + ctData5(j).beamF;

ctData5(j).total15chkv = ctData5(j).ACQ7 + ctData5(j).ACQ1+ ctData5(j).ACQ6;
ctData5(j).total15dose = ctData5(j).beamB + ctData5(j).beamA+ ctData5(j).beamG;


%% PT 57 %%% 
figure
tiledlayout(3,4, 'TileSpacing','tight')
j = 21;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ5); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam E')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ6);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam F')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ1);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam G')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ2);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam H')); set(gca,'FontSize', 18)
%Dose
nexttile
imagesc(medfilt2(ctData5(j).beamE, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam E - 6x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamF, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam F - 15x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamG, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam G - 6x RPO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamH, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam H - 15x RPO')); set(gca,'FontSize', 18)

% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ5); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (medfilt2(ctData5(j).beamE, [30, 30]).* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam E Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ6); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (medfilt2(ctData5(j).beamF, [30, 30]) .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam F Norm')); set(gca,'FontSize', 18)
nexttile
bn3 = double(ctData5(j).ACQ1); msk=double(imbinarize(bn3, 0.2.*prctile(bn3(:), 99.9)));
bn3 = (bn3.*msk) ./ (medfilt2(ctData5(j).beamG, [30, 30]) .* msk);
imagesc(bn3);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam G Norm')); set(gca,'FontSize', 18)
nexttile
bn4 = double(ctData5(j).ACQ2); msk=double(imbinarize(bn4, 0.2.*prctile(bn4(:), 99.9)));
bn4 = (bn4.*msk) ./ (medfilt2(ctData5(j).beamH, [30, 30]) .* msk);
imagesc(bn4);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam H Norm')); set(gca,'FontSize', 18)

ctData5(j).RPO6Xdn = bn3;
ctData5(j).RPO15Xdn = bn4;
ctData5(j).LAO6Xdn = bn1;
ctData5(j).LAO15Xdn = bn2;

j=21;
ctData5(j).total6chkv = ctData5(j).ACQ5 + ctData5(j).ACQ1;
ctData5(j).total6dose = ctData5(j).beamE + ctData5(j).beamG;

ctData5(j).total15chkv = ctData5(j).ACQ6 + ctData5(j).ACQ2;
ctData5(j).total15dose = ctData5(j).beamF + ctData5(j).beamH;

%% PT 61 %%%
figure
tiledlayout(3,2, 'TileSpacing','tight')
j = 24;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ1); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam A')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ6);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam D')); set(gca,'FontSize', 18)
%Dose
nexttile
imagesc(medfilt2(ctData5(j).beamA, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam A - 15x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamD, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam D - 15x RPO')); set(gca,'FontSize', 18)

% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ1); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (medfilt2(ctData5(j).beamA, [30, 30]).* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam A Norm')) ; set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ6); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (medfilt2(ctData5(j).beamD, [30, 30]) .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam D Norm')); set(gca,'FontSize', 18)

ctData5(j).RPO6Xdn = [];
ctData5(j).RPO15Xdn = bn2;
ctData5(j).LAO6Xdn = [];
ctData5(j).LAO15Xdn = bn1 ;
%% PT 65 %%% messed up
figure
tiledlayout(3,3, 'TileSpacing','tight')
j = 25;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ1); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam C')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ6);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam D')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ7);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam D1')); set(gca,'FontSize', 18)
%Dose
nexttile
imagesc(medfilt2(ctData5(j).beamC, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam C - 6x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamD, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam D - 15x RPO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamD1, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam D1 - 6x RPO')); set(gca,'FontSize', 18)

% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ1); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (medfilt2(ctData5(j).beamC, [30, 30]).* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam C Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ6); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (medfilt2(ctData5(j).beamD, [30, 30]) .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam D Norm')); set(gca,'FontSize', 18)
nexttile
bn3 = double(ctData5(j).ACQ7); msk=double(imbinarize(bn3, 0.2.*prctile(bn3(:), 99.9)));
bn3 = (bn3.*msk) ./ (medfilt2(ctData5(j).beamD1, [30, 30]) .* msk);
imagesc(bn3);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam D1 Norm')); set(gca,'FontSize', 18)

ctData5(j).RPO6Xdn = bn3;
ctData5(j).RPO15Xdn = bn2;
ctData5(j).LAO6Xdn = bn1;
ctData5(j).LAO15Xdn = [];

j=25;
ctData5(j).total6chkv = ctData5(j).ACQ1 + ctData5(j).ACQ7;
ctData5(j).total6dose = ctData5(j).beamC + ctData5(j).beamD1;

ctData5(j).total15chkv = ctData5(j).ACQ6;
ctData5(j).total15dose = ctData5(j).beamD;
%% PT 66 %%%
figure
tiledlayout(3,2, 'TileSpacing','tight')
j = 23;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ1); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam E')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ4);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam A')); set(gca,'FontSize', 18)
%Dose
nexttile
imagesc(medfilt2(ctData5(j).beamE, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam E - 6x RPO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamA, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam A - 6x LAO')); set(gca,'FontSize', 18)

% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ1); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (medfilt2(ctData5(j).beamE, [30, 30]).* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam E Norm')) ; set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ4); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (medfilt2(ctData5(j).beamA, [30, 30]) .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam A Norm')); set(gca,'FontSize', 18)

ctData5(j).RPO6Xdn = bn1;
ctData5(j).RPO15Xdn = [];
ctData5(j).LAO6Xdn = bn2;
ctData5(j).LAO15Xdn = [] ;

ctData5(j).total6chkv = ctData5(j).ACQ1 + ctData5(j).ACQ4;
ctData5(j).total6dose = ctData5(j).beamE + ctData5(j).beamA;

%% PT 68 %%% 
figure
tiledlayout(3,5, 'TileSpacing','tight')
j = 18;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ1); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam A')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ2);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam B')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ5);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam E')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ6);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam F')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ7);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam G')); set(gca,'FontSize', 18)
%Dose
nexttile
imagesc(medfilt2(ctData5(j).beamA, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam A - 15x RPO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamB, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam B - 15x RPO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamE, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam E - 6x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamF, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam F - 6x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamG, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam G - 15x LAO')); set(gca,'FontSize', 18)


% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ1); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (medfilt2(ctData5(j).beamA, [30, 30]).* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam A Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ2); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (medfilt2(ctData5(j).beamB, [30, 30]) .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam B Norm')); set(gca,'FontSize', 18)
nexttile
bn3 = double(ctData5(j).ACQ5); msk=double(imbinarize(bn3, 0.2.*prctile(bn3(:), 99.9)));
bn3 = (bn3.*msk) ./ (medfilt2(ctData5(j).beamE, [30, 30]) .* msk);
imagesc(bn3);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam E Norm')); set(gca,'FontSize', 18)
nexttile
bn4 = double(ctData5(j).ACQ6); msk=double(imbinarize(bn4, 0.2.*prctile(bn4(:), 99.9)));
bn4 = (bn4.*msk) ./ (medfilt2(ctData5(j).beamF, [30, 30]) .* msk);
imagesc(bn4);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam F Norm')); set(gca,'FontSize', 18)
nexttile
bn5 = double(ctData5(j).ACQ7); msk=double(imbinarize(bn5, 0.2.*prctile(bn5(:), 99.9)));
bn5 = (bn5.*msk) ./ (medfilt2(ctData5(j).beamG, [30, 30]) .* msk);
imagesc(bn5);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam G Norm')); set(gca,'FontSize', 18)

ctData5(j).RPO6Xdn = [];
ctData5(j).RPO15Xdn = bn1 + bn2;
ctData5(j).LAO6Xdn = bn3 + bn4;
ctData5(j).LAO15Xdn = bn5;


ctData5(j).total6chkv = ctData5(j).ACQ5 + ctData5(j).ACQ6;
ctData5(j).total6dose = ctData5(j).beamE + ctData5(j).beamF;

ctData5(j).total15chkv = ctData5(j).ACQ7 + ctData5(j).ACQ1+ ctData5(j).ACQ2;
ctData5(j).total15dose = ctData5(j).beamB + ctData5(j).beamA+ ctData5(j).beamG;
%% PT 70 %%% 
figure
tiledlayout(3,3, 'TileSpacing','tight')
j = 25;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ1); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam A')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ2);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam B')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ3);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam B1')); set(gca,'FontSize', 18)

%Dose
nexttile
imagesc(medfilt2(ctData5(j).beamA, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam A - 6x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamB, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam B - 6x RPO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamB1, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam B1 - 15x RPO')); set(gca,'FontSize', 18)

% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ1); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (medfilt2(ctData5(j).beamA, [30, 30]).* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam A Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ2); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (medfilt2(ctData5(j).beamB, [30, 30]) .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam B Norm')); set(gca,'FontSize', 18)
nexttile
bn3 = double(ctData5(j).ACQ3); msk=double(imbinarize(bn3, 0.2.*prctile(bn3(:), 99.9)));
bn3 = (bn3.*msk) ./ (medfilt2(ctData5(j).beamB1, [30, 30]) .* msk);
imagesc(bn3);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam B1 Norm')); set(gca,'FontSize', 18)


ctData5(j).RPO6Xdn = bn2;
ctData5(j).RPO15Xdn = bn3;
ctData5(j).LAO6Xdn = bn1;
ctData5(j).LAO15Xdn = [];


ctData5(j).total6chkv = ctData5(j).ACQ1+ctData5(j).ACQ2;
ctData5(j).total6dose = ctData5(j).beamA +ctData5(j).beamB;

ctData5(j).total15chkv = ctData5(j).ACQ3;
ctData5(j).total15dose = ctData5(j).beamB1;
%% PT 72 %%% 
figure
tiledlayout(3,3, 'TileSpacing','tight')
j = 26;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ4); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam A')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ5);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam A1')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ2);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam B')); set(gca,'FontSize', 18)

%Dose
nexttile
imagesc(medfilt2(ctData5(j).beamA, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam A - 6x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamA1, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam A1 - 15x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamB, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam B - 15x RPO')); set(gca,'FontSize', 18)

% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ4); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (medfilt2(ctData5(j).beamA, [30, 30]).* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam A Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ5); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (medfilt2(ctData5(j).beamA1, [30, 30]) .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam A1 Norm')); set(gca,'FontSize', 18)
nexttile
bn3 = double(ctData5(j).ACQ2); msk=double(imbinarize(bn3, 0.2.*prctile(bn3(:), 99.9)));
bn3 = (bn3.*msk) ./ (medfilt2(ctData5(j).beamB, [30, 30]) .* msk);
imagesc(bn3);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam B Norm')); set(gca,'FontSize', 18)


ctData5(j).RPO6Xdn = [];
ctData5(j).RPO15Xdn = bn3;
ctData5(j).LAO6Xdn = bn1;
ctData5(j).LAO15Xdn = bn2;


ctData5(j).total6chkv = ctData5(j).ACQ4;
ctData5(j).total6dose = ctData5(j).beamA;

ctData5(j).total15chkv = ctData5(j).ACQ5+ctData5(j).ACQ2;
ctData5(j).total15dose = ctData5(j).beamA1+ ctData5(j).beamB;
%% PT 76 %%%
figure
tiledlayout(3,2, 'TileSpacing','tight')
j = 21;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ3); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam A')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ1);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam B')); set(gca,'FontSize', 18)
%Dose
nexttile
imagesc(medfilt2(ctData5(j).beamA, [30 30])); colormap(cmap1); axis off; clim([0 5e2])
title(strcat(pt, ' Beam A - 6x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamB, [30 30])); colormap(cmap1); axis off; clim([0 5e2])
title(strcat(pt, ' Beam B - 6x RPO')); set(gca,'FontSize', 18)

% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ3); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (medfilt2(ctData5(j).beamA, [30, 30]).* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam A Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ1); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (medfilt2(ctData5(j).beamB, [30, 30]) .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam B Norm')); set(gca,'FontSize', 18)

ctData5(j).RPO6Xdn = bn2;
% ctData5(j).RPO15Xdn = [];
ctData5(j).LAO6Xdn = bn1;
% ctData5(j).LAO15Xdn = [];


ctData5(j).total6chkv = ctData5(j).ACQ3+ ctData5(j).ACQ1;
ctData5(j).total6dose = ctData5(j).beamA + ctData5(j).beamB;
%% PT 78%%% 
figure
tiledlayout(3,4, 'TileSpacing','tight')
j = 28;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ4); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam A')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ5);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam B')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ2);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam C')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ3);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam D')); set(gca,'FontSize', 18)
%Dose
nexttile
imagesc(medfilt2(ctData5(j).beamA, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam A - 6x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamB, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam B - 15x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamC, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam C - 6x RPO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamD, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam D - 15x RPO')); set(gca,'FontSize', 18)

% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ4); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (medfilt2(ctData5(j).beamA, [30, 30]).* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam A Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ5); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (medfilt2(ctData5(j).beamB, [30, 30]) .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam B Norm')); set(gca,'FontSize', 18)
nexttile
bn3 = double(ctData5(j).ACQ2); msk=double(imbinarize(bn3, 0.2.*prctile(bn3(:), 99.9)));
bn3 = (bn3.*msk) ./ (medfilt2(ctData5(j).beamC, [30, 30]) .* msk);
imagesc(bn3);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam C Norm')); set(gca,'FontSize', 18)
nexttile
bn4 = double(ctData5(j).ACQ3); msk=double(imbinarize(bn4, 0.2.*prctile(bn4(:), 99.9)));
bn4 = (bn4.*msk) ./ (medfilt2(ctData5(j).beamD, [30, 30]) .* msk);
imagesc(bn4);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam D Norm')); set(gca,'FontSize', 18)

ctData5(j).RPO6Xdn = bn3;
ctData5(j).RPO15Xdn = bn4;
ctData5(j).LAO6Xdn = bn1;
ctData5(j).LAO15Xdn = bn2;


ctData5(j).total6chkv = ctData5(j).ACQ4 + ctData5(j).ACQ2;
ctData5(j).total6dose = ctData5(j).beamA + ctData5(j).beamC;

ctData5(j).total15chkv = ctData5(j).ACQ5 + ctData5(j).ACQ3;
ctData5(j).total15dose = ctData5(j).beamB + ctData5(j).beamD;
%% PT 80 %%% 
figure
tiledlayout(3,3, 'TileSpacing','tight')
j = 24;
pt = ctData5(j).Patient;
%Cherenkov
cmax = 3e5;
nexttile
imagesc(ctData5(j).ACQ3); colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam A')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ2);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam B')); set(gca,'FontSize', 18)
nexttile
imagesc(ctData5(j).ACQ1);  colormap(cmap1); axis off; clim([0 cmax])
title(strcat(pt, ' Beam B1')); set(gca,'FontSize', 18)

%Dose
nexttile
imagesc(medfilt2(ctData5(j).beamA, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam A - 6x LAO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamB, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam B - 6x RPO')); set(gca,'FontSize', 18)
nexttile
imagesc(medfilt2(ctData5(j).beamB1, [30 30])); colormap(cmap1); axis off; clim([0 200])
title(strcat(pt, ' Beam B1 - 15x RPO')); set(gca,'FontSize', 18)

% normalize
dnmax = 4e3;
nexttile
bn1 = double(ctData5(j).ACQ3); msk=double(imbinarize(bn1, 0.2.*prctile(bn1(:), 99.9)));
bn1 = (bn1.*msk) ./ (medfilt2(ctData5(j).beamA, [30, 30]).* msk);
imagesc(bn1);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam A Norm')); set(gca,'FontSize', 18)
nexttile
bn2 = double(ctData5(j).ACQ2); msk=double(imbinarize(bn2, 0.2.*prctile(bn2(:), 99.9)));
bn2 = (bn2.*msk) ./ (medfilt2(ctData5(j).beamB, [30, 30]) .* msk);
imagesc(bn2);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam B Norm')); set(gca,'FontSize', 18)
nexttile
bn3 = double(ctData5(j).ACQ1); msk=double(imbinarize(bn3, 0.2.*prctile(bn3(:), 99.9)));
bn3 = (bn3.*msk) ./ (medfilt2(ctData5(j).beamB1, [30, 30]) .* msk);
imagesc(bn3);colormap(cmap1); axis off; clim([0 dnmax])
title(strcat(pt, ' Beam B1 Norm')); set(gca,'FontSize', 18)


ctData5(j).RPO6Xdn = bn2;
ctData5(j).RPO15Xdn = bn3;
ctData5(j).LAO6Xdn = bn1;
ctData5(j).LAO15Xdn = [];


ctData5(j).total6chkv = ctData5(j).ACQ3+ctData5(j).ACQ2;
ctData5(j).total6dose = ctData5(j).beamA+ctData5(j).beamB;

ctData5(j).total15chkv = ctData5(j).ACQ1;
ctData5(j).total15dose = ctData5(j).beamB1;
%% mean CT # 

for i = [2, 5, 6,10,17,19, 21,25]
     figure
     imagesc(ctData5(i).CTnotape); clim([-200 100]); colormap gray; axis off
     h= imellipse(); m1 = createMask(h);
     close
     temp1 = (ctData5(i).CTnotape + 150) .*m1;
    temp1(isnan(temp1))=0;
     ctData5(i).meanCT = mean(nonzeros(temp1(:))) - 150;
end


%% show RPO 15X images



% cints_15rpo = [];
% ct_15rpo = [];          
figure; 
count=1;
for j =[1,3,5,7,12, 17,19,21,25];
    figure;
    imagesc(ctData5(j).CTnotape); clim([-150 150]); colormap gray; axis off
    h= imellipse(); m1 = createMask(h);
    close
    temp1 = (ctData5(j).CTnotape + 150) .*m1;
    temp1(isnan(temp1))=0;
    meanCT =  mean(nonzeros(temp1(:))) - 150;
    figure;
    imagesc(ctData5(j).RPO15Xdn); clim([0 3e3]); colormap(cmap1); axis off
    title(strcat(ctData5(j).Patient)); set(gca,'FontSize', 18)
    h= imellipse(); m2 = createMask(h);
    close
    temp2 = ctData5(j).RPO15Xdn .*m2; ;
    temp2(isnan(temp2))=0;
    cints_15rpo3(count) =  mean(nonzeros(temp2(:))); 
    ct_15rpo3(count) = meanCT;
    hold on
    plot(ct_15rpo3(count),cints_15rpo3(count),'o', 'MarkerSize', 10)

    count = count+1;

end

%  figure; 
% plot(ct_15rpo2, cints_15rpo2,'o', 'MarkerSize', 10, 'MarkerFaceColor', 'b')
% ylabel('Median Chkv/Dose')
% xlabel('Median CT #')
% set(gca, 'FontSize', 20)  
% j ={'1','3','5','7','12','17','19','21','25'};
% legend(j)
%% 
figure;
h1 = plot(rpo15_ct(:,1), 1.0210.*rpo15_chkv(:,1),'o', 'MarkerSize', 10, 'LineWidth',3); hold on;
h2 = plot(rpo15_ct(:,2), 1.0201.*rpo15_chkv(:,2),'o', 'MarkerSize', 10, 'LineWidth',3); hold on;
h3 = plot(rpo15_ct(:,3), 1.0761.*rpo15_chkv(:,3),'o', 'MarkerSize', 10, 'LineWidth',3); hold on;
h4 = plot(rpo15_ct(:,4), 1.0167.*rpo15_chkv(:,4),'o', 'MarkerSize', 10, 'LineWidth',3); hold on;
h5 = plot(rpo15_ct(:,5), 1.0504.*rpo15_chkv(:,5),'o', 'MarkerSize', 10, 'LineWidth',3); hold on;
h6 = plot(rpo15_ct(:,6), 1.0523*rpo15_chkv(:,6),'o', 'MarkerSize', 10, 'LineWidth',3); hold on;
h7 = plot(rpo15_ct(:,7), 1.0658.*rpo15_chkv(:,7),'o', 'MarkerSize', 10, 'LineWidth',3); hold on;
h8 = plot(rpo15_ct(:,8), 1.0648.*rpo15_chkv(:,8),'o', 'MarkerSize', 10, 'LineWidth',3); hold on;
h9 = plot(rpo15_ct(:,9), 1.0761.*rpo15_chkv(:,9),'o', 'MarkerSize', 10, 'LineWidth',3); hold on;

h10 = plot(lao15_ct(:,1), 1.0210.*lao15_chkv(:,1),'diamond', 'MarkerSize', 10, 'LineWidth',3); hold on;
h11 = plot(lao15_ct(:,2), 1.0201.*lao15_chkv(:,2),'diamond', 'MarkerSize', 10, 'LineWidth',3); hold on;
h12 = plot(lao15_ct(:,3), 1.0504.*lao15_chkv(:,3),'diamond', 'MarkerSize', 10, 'LineWidth',3); hold on;
h13 = plot(lao15_ct(:,4), 1.0523*lao15_chkv(:,4),'diamond', 'MarkerSize', 10, 'LineWidth',3); hold on;
h14 = plot(lao15_ct(:,5), 1.0658.*lao15_chkv(:,5),'diamond', 'MarkerSize', 10, 'LineWidth',3); hold on;
h15 = plot(lao15_ct(:,6), 1.0648.*lao15_chkv(:,6),'diamond', 'MarkerSize', 10, 'LineWidth',3); hold on;
ylabel('Mean Chkv/Dose')
xlabel('Mean CT #')
set(gca, 'FontSize', 20)  
legend([h1(1), h2(1), h3(1), h4(1),h5(1),h6(1),h7(1),h8(1),h9(1),h10(1),h11(1),h12(1),h13(1),h14(1),h15(1)],...
    {'PT1 - RPO','PT14 - RPO', 'PT17 - RPO', 'PT3 - RPO', 'PT39 - RPO', 'PT53 - RPO', 'PT55 - RPO', ...
    'PT57 - RPO', 'PT65 - RPO', 'PT1 - LAO','PT14 - LAO', 'PT39 - LAO','PT53 - LAO', 'PT55 - LAO'...
    'PT57 - LAO'}, 'Location', 'eastoutside')

% lao [1,3,12, 17,19,21];
% x = linspace(-200, 100, 300);
% y = ct_fit.p1 .*x + ct_fit.p2;
% hold on 
% plot(x, y, 'b--', 'LineWidth', 2)
% set(gca, 'FontSize', 24)
% xlim([-200 50])
% str=['R^2 =',num2str(ct_goodness.rsquare),newline,...
% 'y = ',num2str(ct_fit.p1),'x+',num2str(ct_fit.p2)]
% t=annotation('textbox',[.15 .9 0 0],'string',str,'FitBoxToText','on', 'FontSize', 18, 'EdgeColor', 'b','LineWidth', 1.5)
% t.FontSize = 20;
%% 
figure;
h1 = scatter(rpo15_ctT, rpo15_chkvT_C,'filled', 'SizeData', 200); hold on;
h2 = scatter(lao15_ctT, lao15_chkvT_C,'filled', 'SizeData', 200); hold on;
alpha(h1, 0.5)
alpha(h2, 0.5)
ylabel('Mean Chkv/Dose')
xlabel('Mean CT #')

% 
x = linspace(-150, 0, 300);
y = ct15C_fit.p1 .*x + ct15C_fit.p2;
hold on 
plot(x, y, 'k--', 'LineWidth', 2)
set(gca, 'FontSize', 24)
xlim([-130 0])
ylim([0 2500])
str=['R^2 =',num2str(ct15C_goodness.rsquare),newline,...
'y = ',num2str(ct15C_fit.p1),'x+',num2str(ct15C_fit.p2)]
t=annotation('textbox',[.15 0.3 0 0],'string',str,'FitBoxToText','on', 'FontSize', 18, 'EdgeColor', 'k','LineWidth', 1)
t.FontSize = 20;
set(gca, 'FontSize', 20)  
legend({'RPO' ,'LAO', 'fit'})



%% show LAO 15X images
  %j =[1,3,4,6,12,17,19,21];


         
figure; 
count=1;
for j =[1,3,12, 17,19,21];
    figure;
    imagesc(ctData5(j).CTnotape); clim([-150 150]); colormap gray; axis off
    h= imellipse(); m1 = createMask(h);
    close
    temp1 = (ctData5(j).CTnotape + 150) .*m1;
    temp1(isnan(temp1))=0;
    meanCT =  mean(nonzeros(temp1(:))) - 150;
    figure;
    imagesc(ctData5(j).LAO15Xdn); clim([0 3e3]); colormap(cmap1); axis off
    title(strcat(ctData5(j).Patient)); set(gca,'FontSize', 18)
    h= imellipse(); m2 = createMask(h);
    close
    temp2 = ctData5(j).LAO15Xdn .*m2; ;
    temp2(isnan(temp2))=0;
    cints_15lao1(count) =  mean(nonzeros(temp2(:))); 
    ct_15lao1(count) = meanCT;
    hold on
    plot(ct_15lao1(count),cints_15lao1(count),'o', 'MarkerSize', 10)

    count = count+1;

end

%  figure; 
% plot(ct_15rpo2, cints_15rpo2,'o', 'MarkerSize', 10, 'MarkerFaceColor', 'b')
% ylabel('Median Chkv/Dose')
% xlabel('Median CT #')
% set(gca, 'FontSize', 20)  
% j ={'1','3','5','7','12','17','19','21','25'};
% legend(j)
%% show RPO 6X images

      
figure; 
count=1;
for j =[1,10,17,19,21,25];
    figure;
    imagesc(ctData5(j).CTnotape); clim([-150 150]); colormap gray; axis off
    h= imellipse(); m1 = createMask(h);
    close
    temp1 = (ctData5(j).CTnotape + 150) .*m1;
    temp1(isnan(temp1))=0;
    meanCT =  mean(nonzeros(temp1(:))) - 150;
    figure;
    imagesc(ctData5(j).RPO6Xdn); clim([0 3e3]); colormap(cmap1); axis off
    title(strcat(ctData5(j).Patient)); set(gca,'FontSize', 18)
    h= imellipse(); m2 = createMask(h);
    close
    temp2 = ctData5(j).RPO6Xdn .*m2; ;
    temp2(isnan(temp2))=0;
    cints_6rpo1(count) =  mean(nonzeros(temp2(:))); 
    ct_6rpo1(count) = meanCT;
    hold on
    plot(ct_6rpo1(count),cints_6rpo1(count),'o', 'MarkerSize', 10)

    count = count+1;

end

%  figure; 
% plot(ct_15rpo2, cints_15rpo2,'o', 'MarkerSize', 10, 'MarkerFaceColor', 'b')
% ylabel('Median Chkv/Dose')
% xlabel('Median CT #')
% set(gca, 'FontSize', 20)  
% j ={'1','3','5','7','12','17','19','21','25'};
% legend(j)

%% show LAO 6X images
  %j =[1,3,4,6,12,17,19,21];


         

count=1;
for j =[1,3,5,7,10,12, 17,19,21,25];
    figure;
    imagesc(ctData5(j).CTnotape); clim([-150 150]); colormap gray; axis off
    h= imellipse(); m1 = createMask(h);
    close
    temp1 = (ctData5(j).CTnotape + 150) .*m1;
    temp1(isnan(temp1))=0;
    meanCT =  mean(nonzeros(temp1(:))) - 150;
    figure;
    imagesc(ctData5(j).LAO6Xdn); clim([0 3e3]); colormap(cmap1); axis off
    title(strcat(ctData5(j).Patient)); set(gca,'FontSize', 18)
    h= imellipse(); m2 = createMask(h);
    close
    temp2 = ctData5(j).LAO6Xdn .*m2; ;
    temp2(isnan(temp2))=0;
    cints_6lao3(count) =  mean(nonzeros(temp2(:))); 
    ct_6lao3(count) = meanCT;
    hold on
    plot(ct_6lao3(count),cints_6lao3(count),'o', 'MarkerSize', 10)

    count = count+1;

end

%  figure; 
% plot(ct_15rpo2, cints_15rpo2,'o', 'MarkerSize', 10, 'MarkerFaceColor', 'b')
% ylabel('Median Chkv/Dose')
% xlabel('Median CT #')
% set(gca, 'FontSize', 20)  
% j ={'1','3','5','7','12','17','19','21','25'};
% legend(j)

%% 
figure;
h1 = plot(rpo6_ct(:,1), rpo6_chkv(:,1),'o', 'MarkerSize', 10, 'LineWidth',3); hold on;
h2 = plot(rpo6_ct(:,2), rpo6_chkv(:,2),'o', 'MarkerSize', 10, 'LineWidth',3); hold on;
h3 = plot(rpo6_ct(:,3), rpo6_chkv(:,3),'o', 'MarkerSize', 10, 'LineWidth',3); hold on;
h4 = plot(rpo6_ct(:,4), rpo6_chkv(:,4),'o', 'MarkerSize', 10, 'LineWidth',3); hold on;
h5 = plot(rpo6_ct(:,5), rpo6_chkv(:,5),'o', 'MarkerSize', 10, 'LineWidth',3); hold on;
h6 = plot(rpo6_ct(:,6), rpo6_chkv(:,6),'o', 'MarkerSize', 10, 'LineWidth',3); hold on;

h7 = plot(lao6_ct(:,1), lao6_chkv(:,1),'diamond', 'MarkerSize', 10, 'LineWidth',3); hold on;
h8 = plot(lao6_ct(:,2), lao6_chkv(:,2),'diamond', 'MarkerSize', 10, 'LineWidth',3); hold on;
h9 = plot(lao6_ct(:,3), lao6_chkv(:,3),'diamond', 'MarkerSize', 10, 'LineWidth',3); hold on;
h10 = plot(lao6_ct(:,4), lao6_chkv(:,4),'diamond', 'MarkerSize', 10, 'LineWidth',3); hold on;
h11 = plot(lao6_ct(:,5), lao6_chkv(:,5),'diamond', 'MarkerSize', 10, 'LineWidth',3); hold on;
h12 = plot(lao6_ct(:,6), lao6_chkv(:,6),'diamond', 'MarkerSize', 10, 'LineWidth',3); hold on;
h13 = plot(lao6_ct(:,7), lao6_chkv(:,7),'diamond', 'MarkerSize', 10, 'LineWidth',3); hold on;
h14 = plot(lao6_ct(:,8), lao6_chkv(:,8),'diamond', 'MarkerSize', 10, 'LineWidth',3); hold on;
h15 = plot(lao6_ct(:,9), lao6_chkv(:,9),'diamond', 'MarkerSize', 10, 'LineWidth',3); hold on;
h16 = plot(lao6_ct(:,10), lao6_chkv(:,10),'diamond', 'MarkerSize', 10, 'LineWidth',3); hold on;
ylabel('Median Chkv/Dose')
xlabel('Median CT #')
set(gca, 'FontSize', 20)  
legend([h1(1), h2(1), h3(1), h4(1),h5(1),h6(1),h7(1),h8(1),h9(1),h10(1),h11(1),h12(1),h13(1),h14(1),h15(1), h16(1)],...
    {'PT1 - RPO','PT35 - RPO', 'PT53 - RPO', 'PT55 - RPO', 'PT57 - RPO', 'PT65 - RPO', 'PT1 - LAO', ...
    'PT14 - LAO', 'PT17 - LAO', 'PT3 - LAO','PT35 - LAO', 'PT39 - LAO','PT53 - LAO', 'PT55 - LAO'...
    'PT57 - LAO','PT65 - LAO'})
% 
% x = linspace(-200, 100, 300);
% y = ct_fit.p1 .*x + ct_fit.p2;
% hold on 
% plot(x, y, 'b--', 'LineWidth', 2)
% set(gca, 'FontSize', 24)
% xlim([-200 50])
% str=['R^2 =',num2str(ct_goodness.rsquare),newline,...
% 'y = ',num2str(ct_fit.p1),'x+',num2str(ct_fit.p2)]
% t=annotation('textbox',[.15 .9 0 0],'string',str,'FitBoxToText','on', 'FontSize', 18, 'EdgeColor', 'b','LineWidth', 1.5)
% t.FontSize = 20;
%% 
figure;
h1 = scatter(rpo6_ctT, rpo6_chkvT_C,'filled', 'SizeData', 200); hold on;
h2 = scatter(lao6_ctT, lao6_chkvT_C,'filled', 'SizeData', 200); hold on;
alpha(h1, 0.5)
alpha(h2, 0.5)
ylabel('Median Chkv/Dose')
xlabel('Median CT #')

% 
x = linspace(-150, 0, 300);
y = ct6C_fit.p1 .*x + ct6C_fit.p2;
hold on 
plot(x, y, 'k--', 'LineWidth', 2)
set(gca, 'FontSize', 24)
xlim([-150 0])
str=['R^2 =',num2str(ct6C_goodness.rsquare),newline,...
'y = ',num2str(ct6C_fit.p1),'x+',num2str(ct6C_fit.p2)]
t=annotation('textbox',[.15 0.3 0 0],'string',str,'FitBoxToText','on', 'FontSize', 18, 'EdgeColor', 'k','LineWidth', 1)
t.FontSize = 20;
set(gca, 'FontSize', 20)  
legend({'RPO' ,'LAO', 'fit'})


