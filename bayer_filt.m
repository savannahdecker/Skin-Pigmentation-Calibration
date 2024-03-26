% Bayer spatial median filter (before demosaicing)

function filtered_bayer = bayer_filt(im, filt_size)

%im is a bayer filter image with an rggb pattern and filt_size is a number
%that will size a square median filter of that size 

   
    bayerPattern = [1 2; 2 3];
    
    % for rggb pattern
    redChannel = im(1:2:end, 1:2:end);
    greenChannel1 = im(1:2:end, 2:2:end);
    greenChannel2 = im(2:2:end, 1:2:end);
    blueChannel = im(2:2:end, 2:2:end);
    
    filt = [filt_size, filt_size];
    
    filteredRed = medfilt2(redChannel, filt);
    filteredGreen1 = medfilt2(greenChannel1, filt);
    filteredGreen2 = medfilt2(greenChannel2, filt);
    filteredBlue = medfilt2(blueChannel, filt);


    % Combine the filtered color channels into a single image
    filtered_bayer = zeros(size(im));
    filtered_bayer(1:2:end, 1:2:end) = filteredRed;
    filtered_bayer(1:2:end, 2:2:end) = filteredGreen1;
    filtered_bayer(2:2:end, 1:2:end) = filteredGreen2;
    filtered_bayer(2:2:end, 2:2:end) = filteredBlue;
end

