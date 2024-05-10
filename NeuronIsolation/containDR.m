function drFlag = containDR(imgNeuron_store)
% This function checks if the neurons contain darkRegion
% Input 16-bit neuron image and return 1 if there is a relatively large
% dark region in the image, otherwise return 0 
% Minimum detectable dark region size is 1% of the original image
    darkRegionBW = double(imgNeuron_store)==0;
    darkRegionDimX = [];
    darkRegionDimY = [];
    CC_dark = bwconncomp(darkRegionBW);
    stats = regionprops(CC_dark,'BoundingBox'); % Stats format: [ x0 y0(upper-left corner) x-width y-width] 
    imgSize = CC_dark.ImageSize;
    for i = 1:CC_dark.NumObjects
        darkRegionDimX = [darkRegionDimX;  stats(i).BoundingBox(3)];
        darkRegionDimY = [darkRegionDimY;  stats(i).BoundingBox(4)];
    end
    if any(darkRegionDimX >= imgSize(2)/10) || any(darkRegionDimY >= imgSize(1)/10)
        drFlag = 1;
    else 
        drFlag = 0;
    end
end