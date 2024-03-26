function output = roi_mean(ims, varargin)
  
%this function will output the mean values of the ROI selected, for a
%number of ROIs, n

%ims is a cell containing all images where the same ROI should be analyzed

% varargin could be the number of regions in a single image
% you want to analyze. if excluded, n is assumed to be 1

addpath('/Users/f004mg7/Library/CloudStorage/GoogleDrive-savannah.m.decker.th@dartmouth.edu/My Drive/Research/Scripts')
load cmap1.mat cmap1

    if nargin >1
        n = varargin{1};
    else 
        n =1;
    end
    
    [r, c] = size(ims);
    
 
% define the same roi for all images
    if n == 1
        output = [];
        figure; imagesc(ims{1,1}); cmax = 1.5.*prctile(ims{1,1}(:), 99);
        clim([0 cmax]); 
        colormap(cmap1)
        h = imrect(); m = createMask(h);
        % apply the same mask to each image in ims
        for k = 1:c
            im = ims{1,k};
            temp = double(im).*double(m);
            output_temp = mean(nonzeros(temp(:)), 'omitnan');
            output(k) = output_temp;
        end
        close
 % nested for loop if dealing with multiple ROIs in multiple images
    else 
        output = [];
        hs = {}; ms = {};
        figure; imagesc(ims{1,1}); cmax = 1.5.*prctile(ims{1,1}(:), 99);
        clim([0 cmax]); colormap(cmap1)
        for i = 1:n
            hs{i} = imrect(); ms{i} = createMask(hs{i});
            for j = 1:c
                im = ims{1, j};
                temp = double(im).*double(ms{i});
                [g, f, z] = size(temp);
                for  k = 1:z
                    temp2 = temp(:,:,k);
                    output_temp = mean(nonzeros(temp2(:)), 'omitnan');
                    output(i, j, k) = output_temp;
                end
            end
            %close 
        end
        close 
         
    end
    
end
