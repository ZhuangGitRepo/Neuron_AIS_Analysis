function loopFlag = containLoop(neuronBW2,loopPeriThresh)
% This function checks if the neurons contain loops
% Input binary neuron image and return the flag 0 if no loop, 1 if there is
% a loop (minimum detectable loop perimeter is specified by loopPeriThresh)
loopRegion =  imclearborder(~neuronBW2);
loopPerimeter = regionprops(loopRegion ,'Perimeter');
loopPerimeter = cat(1,loopPerimeter.Perimeter);
loopFlag = any(loopPerimeter>loopPeriThresh);
end

