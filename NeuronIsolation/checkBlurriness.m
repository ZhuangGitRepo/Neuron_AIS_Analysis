function [blurriness blurFlag] =  checkBlurriness(imgNeuron_store,blurThresh)
% % This function checks if the input neuron image is blurry
% % Input 16-bit neuron image and return 1 if the blurriness of the image(calculated by 
% % taking fft of the input image and taking the average of the 45% ((2/3)^2) highest frequencies).
% % If the exceed the given threshold, otherwise return 0

[ny nx] = size(imgNeuron_store);

fftData = fft2(imgNeuron_store);
fftDataShift = abs(fftshift(fftData));
Fthresh = max(fftDataShift(:))/1000; % choose one thousandth of the maximum value in the frequency domain.
blurrinessMatrix = fftDataShift > Fthresh;
blurriness = sum(blurrinessMatrix(:))/(ny*nx);

if blurriness < blurThresh
    blurFlag = 1;
else
    blurFlag = 0;
end
end

