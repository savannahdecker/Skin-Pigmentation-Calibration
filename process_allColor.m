% loop through folders, each named after the patient # and and contains 1
% folder of the acquisition with the desired gantry position for color 
% Savannah Decker July 2023

cd('/Users/f004mg7/Downloads/Moffitt/Color Image Data/Right Breast')
addpath('/Users/f004mg7/Library/CloudStorage/GoogleDrive-savannah.m.decker.th@dartmouth.edu/My Drive/Research/Scripts')

% READ IN ALL COLOR IMAGES ANS SAVE IN STRUCT
files = dir(pwd);
tempfiles = {};

for f = 1:size(files,1)
    if startsWith(files(f).name, 'PT')
        tempfiles = [tempfiles;files(f).name]
    end
end

fld = files(1).folder
% color_data = struct;
% 
for i = 1:numel(tempfiles)
    cd(strcat(fld, '/', tempfiles{i,1}))
    temp_fld = dir(pwd);
    cd(temp_fld(end).name)
    temp = read_dovi_c('meas_cam0.dovi', 's1');
    color_data(32).pt = tempfiles{i,1};
    color_data(32).colorIm = temp;
    cd(fld)
end

%% DEMOSAIC AND COLOR CORRECT EACH IMAGE AND SAVE
tic
for i = [32]
    pt_test = color_data(i).colorIm;
    test_NOF = 1;
    filtered_bayer = bayer_filt(pt_test, 3);
    bkg_im =filtered_bayer(1:end-mod(end,2),1:end-mod(end,2));
    bkg_image_de = demosaic(uint32(bkg_im), 'rggb');
    bkg_image_filt(:, :, 1) = medfilt2(bkg_image_de(:, :, 1), [5 5]);%.* cRGB(1);
    bkg_image_filt(:, :, 2) = medfilt2(bkg_image_de(:, :, 2), [5 5]);%.* cRGB(2);
    bkg_image_filt(:, :, 3) = medfilt2(bkg_image_de(:, :, 3), [5 5]);%.* cRGB(3);
    %
    %
    rescale_factor =2;
    color_bkg= rescale(double(bkg_image_filt),0, rescale_factor);
    bkg_cc = imapplymatrix(ccm(1:3,:)',color_bkg,ccm(4,:));
   %figure; imagesc(bkg_cc)
    shadow_lab = rgb2lab(bkg_cc);
    %The values of luminosity span a range from 0 to 100. Scale the values to the range [0 1], which is the expected range of images with data type double.
    
    max_luminosity = 100;
    L = shadow_lab(:,:,1)/max_luminosity;
    %Perform the three types of contrast adjustment on the luminosity channel, and keep the a* and b* channels unchanged. Convert the images back to the RGB color space.
     
    shadow_imadjust = shadow_lab;
    shadow_imadjust(:,:,1) = imadjust(L)*max_luminosity;
    shadow_imadjust = lab2rgb(shadow_imadjust);



    color_data(i).cc_colorIm = shadow_imadjust ;
   % figure;
    % imagesc(shadow_imadjust ); axis off; title(strcat('Color-Corrected  ', color_data(i).pt))
    % set(gca, 'FontSize', 20)
    cd('/Users/f004mg7/Downloads/Moffitt/Color Images/September2023')
   % saveas(gcf, strcat(color_data(i).pt, '.svg'));
   % saveas(gcf, strcat(color_data(i).pt, '.png'));
    imwrite(shadow_imadjust , strcat(color_data(i).pt, '.png'));
    color_data(i).colorPNG = imread(strcat(color_data(i).pt, '.png'));
end
%save('colordata.mat', "color_data")
toc
%
cd('/Users/f004mg7/Downloads/Moffitt/Color Images/September2023')
files = dir(pwd);
tempfiles = {};

for f = 1:size(files,1)
    if startsWith(files(f).name, 'PT')
        tempfiles = [tempfiles;files(f).name]
    end
end

fld = files(1).folder
color_data = struct;

for i = 1:length(tempfiles)
    color_data(i).pt = tempfiles{i,1}(1:end-4);
    color_data(i).colorIm = imread(tempfiles{i,1});
end

figure;
tiledlayout(6,6, "TileSpacing","tight")
for i = 1:length(color_data)
    nexttile
    imagesc(color_data(i).colorIm);
  % colormap(gray)
  %  colorbar; caxis([0 4e3])
    axis off;
    title(strcat(color_data(i).pt));%, RB_pt_data(i).Ethnicty));
end

%%
% for i = 1:length(ctData4)
%     ctData4(i).colorIm = color_data(i).colorIm;
% end
%% check 
% figure;
% tiledlayout(4,6, "TileSpacing","tight")
% for i = 1:length(ctData5)
%     nexttile
%     imagesc(ctData5(i).ACQ1);
% %   colormap(gray); caxis([-200 100])
%     colormap(cmap1); colorbar; caxis([0 2e5])
%     axis off;
%     title(strcat(ctData5(i).Patient));%, RB_pt_data(i).Ethnicty));
% end