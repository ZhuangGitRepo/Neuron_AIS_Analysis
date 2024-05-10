function [somaRadiusList areaList DAPIlist region2D circularityList blurrinessList contactingNeuron unidentifiedObj roundCell loopNum blurImage darkRegion] = detectNeuron(multiMarkersIMG,intThresh,pixelThreshL,...
    pixelThreshH,extendLenth,nucleusPixelThreshL,nucleusPixelThreshH,newfolderPath,neuronLabel,tileLabel,MarkerStrings)

% multiMarkersIMG should be 16 bit, and 1st and 2nd channel of multiMarkersIMG 
% should always be DAPI and tubulin respectively.


% Convert DAPI and tubulin images to 8-bit, since their intensities are
% high (better for APP2 automated tracing).
% For Tau and AnkG, Keep original 16-bit, as their intensities
% are low, coverting to 8-bit image would lose these low level contrast
index_DAPI = find(MarkerStrings == 'DAPI');
index_tubulin = find(MarkerStrings == 'tubulin');
index_3 = find(MarkerStrings == 'AnkG');
index_4 = find(MarkerStrings == 'Tau');

tile_DAPI_8bit =  uint8(multiMarkersIMG(:,:,index_DAPI)/256);
tile_tubulin_8bit =  uint8(multiMarkersIMG(:,:,index_tubulin)/256);
tile_tubulin_16bit = multiMarkersIMG(:,:,index_tubulin);
tile_Marker3_16bit = multiMarkersIMG(:,:,index_3);
tile_Marker4_16bit = multiMarkersIMG(:,:,index_4);


tile_DAPI_8bit = medfilt2(tile_DAPI_8bit);
tile_tubulin_8bit = medfilt2(tile_tubulin_8bit);
tile_Marker3_16bit = medfilt2(tile_Marker3_16bit);
tile_Marker4_16bit = medfilt2(tile_Marker4_16bit);
% tile_tubulin_16bit = medfilt2(tile_tubulin_16bit); % for blurriness
% detection, do not apply this filter to tubulin16bit?

% imgTubulin = tile_tubulin_8bit;
imgTubulin = medfilt2(tile_tubulin_16bit);
imgDAPI = tile_DAPI_8bit;
imgTubulin16bit = tile_tubulin_16bit;
imgMarker3 =  tile_Marker3_16bit;
imgMarker4 = tile_Marker4_16bit;

% detect neurons and extract the range of x- & y- coordinates each neuron covers
somaRadiusList = [];
DAPIlist = [];
areaList = [];

region2D = [];
circularityList = [];
blurrinessList = [];

contactingNeuron = 0;
unidentifiedObj = 0;
roundCell = 0; %  % Number of non-contacting cell with round shape
loopNum = 0; % Number of non-contacting neurons with a loop
blurImage = 0;
darkRegion = 0;



%median filter used to enhance image quality, reduce noice and while
%keeping the boundary sharp,  default setting uses a 3 by 3 window
% imgTubulin = medfilt2(imgTubulin); % 8 bit tubulin image is used for connected components identification
imgDAPI = medfilt2(imgDAPI);
imgMarker3 = medfilt2(imgMarker3);
imgMarker4 = medfilt2(imgMarker4);

% intThresh =  quantile(imgTubulin(:),0.9);
imgTubulinBW2 = imgTubulin>intThresh; %intThresh = 5 too large for data 0758
figure(100)
imshow(imgTubulinBW2)

% Remove small components less than pixelThresh pixels(convert to square microns later)
% imgTubulinBW2 = bwareaopen(imgTubulinBW2,pixelThreshL);
imgTubulinBW2 = bwareafilt(imgTubulinBW2,[pixelThreshL pixelThreshH]); 
% Remove neurons that touch the boundaries of the image, as the morphology is probaby incomplete
imgTubulinBW2 =  imclearborder(imgTubulinBW2);

%Extract connected component from the binary image
CC = bwconncomp(imgTubulinBW2);
% CC_Cnn = shapeCnnBW2;
% CC_Cnn(CC.PixelIdxList{1}) = 0;
% imshow(shapeCnnBW2.*~CC_Cnn)

[imgEdgeY imgEdgeX]  = size(imgTubulin);%detect

for i = 2:CC.NumObjects
    idxCC = cell2mat(CC.PixelIdxList(i));
    [row,col] = ind2sub(size(imgTubulinBW2),idxCC);
    Width = [min(row)-extendLenth max(row)+extendLenth]; %y-range
    Length = [min(col)-extendLenth max(col)+extendLenth]; %x-range
    Width(Width<=0) = 1;
    Width(Width>imgEdgeY) = imgEdgeY;
    Length(Length<=0) = 1;
    Length(Length>imgEdgeX) = imgEdgeX;        
    imgNeuron = imgTubulin(Width(1):Width(2),Length(1):Length(2)); % for plotting 8bit image
    imgNeuron_store = imgTubulin16bit(Width(1):Width(2),Length(1):Length(2)); %16bit scale to get 8bit
%     imgNeuron_store = uint8(double(imgNeuron_store)/double(max(max(imgNeuron_store)))*256);
    imgDAPI_store = imgDAPI(Width(1):Width(2),Length(1):Length(2));
    img3_store = imgMarker3(Width(1):Width(2),Length(1):Length(2));
    img4_store = imgMarker4(Width(1):Width(2),Length(1):Length(2));
    neuronBW2 = false(size(imgTubulin));
%     for i = 1:1:length(row)
%         neuronBW2(row(i),col(i)) = true;
%     end
    
    neuronBW2 = neuronBW2(:);
    neuronBW2(idxCC) = true;
    neuronBW2 = reshape(neuronBW2,size(imgTubulin));
%     figure
%     imshow(neuronBW2)
    neuronBW2 = neuronBW2(Width(1):Width(2),Length(1):Length(2));
    
    
%     figure
%     imshow(neuronBW2)  
    D = bwdist(~neuronBW2);
    if max(max(D))<15 || max(max(D))> 60 %  soma size pixel range 
        display('Warning: soma size too small or too large')
        unidentifiedObj = unidentifiedObj + 1; 
        continue
    end
    r = 7; % radius of disk-shaped dilation structuring element
    SE = strel('disk',r);
    neuronBW2Dilate = imdilate(neuronBW2,SE);
    
%     neuronExtract = imgNeuron_store;
%     imgNeuron_store2 = uint8(double(imgNeuron_store)/double(max(max(imgNeuron_store)))*256);
%     figure(3)
%     imshow(imgNeuron_store2)
%     se = strel('disk',35);
%     background = imopen(imgNeuron_store2,se);
%     figure(4)
% %     imshow(uint8(double(background)/double(max(max(background)))*256))
%     imshow(background)
%     imgNeuron_store3 = imgNeuron_store2 - background;
%     figure(5)
%     imshow(imgNeuron_store3)
%     imgNeuron_store4 = imadjust(imgNeuron_store3,stretchlim(imgNeuron_store3,[0.01 0.998])); %saturating 1% of the data at both low and high intensities and by stretching the intensity values to fill the uint8 dynamic range.
% %     imgNeuron_store4 = imadjust(imgNeuron_store3,[0.01,0.7]);
%     figure(6)
%     imshow(imgNeuron_store4)
% %     imwrite(imgNeuron_store4,'C:\Users\ericx\Documents\testtemp\test11.tiff')
%     idxCC_background = find(neuronBW2Dilate==0);
% %     neuronExtract(idxCC_background) = 0;
%     figure(7)
%     imgNeuron_store5 =  imgNeuron_store4;
%     imgNeuron_store5(idxCC_background) = 0;
%     imshow(imgNeuron_store5)
% %     imwrite(imgNeuron_store5,'C:\Users\ericx\Documents\testtemp\test12.tiff')

    % Image preprocessing to enhance contrast & remove other neurons in the
    % image
    neuronExtract = imgNeuron_store;
    neuronExtract = uint8(double(neuronExtract)/double(max(max(neuronExtract)))*256);
    %se = strel('disk',35); %previous setting
    se = strel('disk',35); % changed on 13/04/21
    background = imopen(neuronExtract,se);
    neuronExtract = neuronExtract - background;
    neuronExtract = imadjust(neuronExtract,stretchlim(neuronExtract,[0.01 0.998])); %saturating 1% of the data at both low and high intensities and by stretching the intensity values to fill the uint8 dynamic range.
    idxCC_background = find(neuronBW2Dilate==0);
    neuronExtract(idxCC_background) = 0;
    
    % Check contacting soma here!!!
    idxCC_Dilate = find(neuronBW2Dilate==1);
    imgNucleus = imgDAPI_store;
    %detect soma by its physical size, diameter of nucleaus should be at least
    %1 or 2 or 3 microns
    DAPIIntThresh = round(mean2(imgNucleus))+2.5; %works for 86  
%     DAPIIntThresh = round(mean2(imgNucleus))+ 5; % 681
%     DAPIIntThresh = round(mean2(imgNucleus))+5; % 0758
    DAPIBW2 = imgNucleus > DAPIIntThresh ;
    %remove small components less than nucleusPixelThreshL (Lower bound)
    DAPIBW2 = bwareaopen(DAPIBW2,nucleusPixelThreshL);
    % imgTubulinBW2 = bwareafilt(imgTubulinBW2,[pixelThreshL pixelThreshH]); 
    

    figure(1)
    subplot(2,3,1)
    imshow(imgNeuron,[])
    title('Neuron images')
    subplot(2,3,2)
    imshow(neuronBW2)
    title('BW Neuron')
    subplot(2,3,3)
    imshow(neuronBW2Dilate)
    title('BW neuron dilation')
    subplot(2,3,4)
    imshow(imgNucleus,[])
    title('DAPI Nucleus')
    subplot(2,3,5)
    imshow(DAPIBW2)
    title('BW Nucleus')
    subplot(2,3,6)
    imshow(neuronExtract,[])
    title('Extracted neuron') 
    
    se = strel('disk',30);
    figure(111111)
    imshow(imtophat(neuronExtract,se),[])
    
    figure(2)
    nuclei_neuronBW = 0.2*neuronBW2Dilate + 0.8*DAPIBW2;
    imshow(nuclei_neuronBW);
    colormap([0 0 0;1 1 0; 1 0 0; 0 0 1])
    nuclei_neuronBW = 1+1*neuronBW2Dilate + 2*DAPIBW2; % Matlab color different integer with different color
    CMap = [0 0 0;1 1 0; 1 0 0; 0 0 1];
    nuclei_neuronBW_RGB = ind2rgb(nuclei_neuronBW,CMap);
    imshow(nuclei_neuronBW_RGB)
    title('Nuclei & Neuron overlap')   
    
    CC_DAPI = bwconncomp(DAPIBW2); % Connected components;
    % detect number of soma in the extracted neuron image, if somaNum is greater than 1, the image should be skipped
    somaNum = 0; 
    for k = 1:1:CC_DAPI.NumObjects
        idxCC_DAPI = cell2mat(CC_DAPI.PixelIdxList(k));
        if sum(ismember(idxCC_DAPI,idxCC_Dilate))>= 1
            somaNum = somaNum + 1;
        end
    end
    
    if somaNum > 1 
        display('Warning: contacting neurons!')
        contactingNeuron = contactingNeuron + 1;
        continue
    elseif somaNum == 0
        display('Warning: No nucleus detected!')
        unidentifiedObj = unidentifiedObj + 1; 
        continue
    elseif sum(DAPIBW2(idxCC_Dilate)) > nucleusPixelThreshH % If physical diameter of nucleus > sqrt(2000/pi)*0.325*2 = 14.5 micron
        display('Warning: Nucleus too large! (Clustered soma)')
        unidentifiedObj = unidentifiedObj + 1; 
        continue
    end
    
    loopFlag  = 0;
    roundFlag = 0;
    drFlag = 0;
    blurFlag = 0;
    
    % Check if image contain loops using neuronBW2 as input 
    % (minimum detectable loop perimeter is specified by loopPeriThresh)
    loopPeriThresh = 60;
    loopFlag = containLoop(neuronBW2,loopPeriThresh);
    if loopFlag == 1
        display('The neuron contans a large loop')
        loopNum = loopNum + 1;
        continue
    end
    
%     % Check if the neuron deteceted is a round object using neuronBW2 as input 
%     % (minimum circularity threshold is specified by circularityThresh)
%     circularityThresh = 0.25; % circularity higher than this threshold is classidied as round
%     [circularity roundFlag] = isRound(neuronBW2,circularityThresh);
%     if roundFlag == 1
%         fprintf('Circularity:%.3f, The identified object is a round cell!\n',circularity)
%         roundCell = roundCell + 1;
%         continue
%     end
    
    % Check if the neuron image contains a relatively large dark region using
    % 16-bit neuron image(tubulin) as input. (Minimum detectable dark region size is 1% of the original image)
    drFlag = containDR(imgNeuron_store);
    if drFlag == 1
        display('Warning: image contain dark regions!')
        darkRegion = darkRegion + 1;
        continue
    end
    
    % Check if the neuron image is blurry using 16-bit neuron image(tubulin) as input
    % (minimum circularity threshold is specified by circularityThresh)
    blurThresh = 0.02; 
    [blurriness blurFlag] =  checkBlurriness(imgNeuron_store,blurThresh);
    if blurFlag == 1
        display('Warning: the image is blurry!')
        blurImage = blurImage + 1;
        continue
    end
      
    neuronLabel = neuronLabel + 1;
%     imgDAPI_store(idxCC_remove) = 0;
%     imgMAP_store(idxCC_remove) = 0;
%     imgAnkG_store(idxCC_remove) = 0;
    pathstr_DAPI = sprintf('%s\\tile%d_DAPI%d.tiff',newfolderPath,tileLabel,neuronLabel);
    pathstr_tubulin = sprintf('%s\\tile%d_tubulin%d.tiff',newfolderPath,tileLabel,neuronLabel);
%         pathstr_MAP = sprintf('%s\\tile%d_MAP%d.tiff',newfolderPath,tileLabel,neuronLabel);
    pathstr_3 = sprintf('%s\\tile%d_%s%d.tiff',newfolderPath,tileLabel,MarkerStrings(index_3),neuronLabel);
    pathstr_4 = sprintf('%s\\tile%d_%s%d.tiff',newfolderPath,tileLabel,MarkerStrings(index_4),neuronLabel);
    pathstr_nuclei_neuron = sprintf('%s\\tile%d_neuron_neuclei_overlap%d.tiff',newfolderPath,tileLabel,neuronLabel);
    
    imwrite(neuronExtract,pathstr_tubulin)
    imwrite(imgDAPI_store,pathstr_DAPI)
    imwrite(img3_store,pathstr_3)
    imwrite(img4_store,pathstr_4)
    imwrite(nuclei_neuronBW_RGB,pathstr_nuclei_neuron)
    
    pathstr_tubOriginal = sprintf('%s\\tile%d_tubOrignal%d.tiff',newfolderPath,tileLabel,neuronLabel);
    imgNeuron_store = medfilt2(imgNeuron_store);
    imwrite(imgNeuron_store,pathstr_tubOriginal)
    
    region2D = [region2D; Width Length];
    somaRadiusList = [somaRadiusList; max(max(D))];
    areaList = [areaList; sum(sum(neuronBW2))];
    DAPIlist = [DAPIlist; sum(DAPIBW2(idxCC_Dilate))]; 
    circularityList = [circularityList; circularity];
    blurrinessList = [blurrinessList; blurriness];
    
    % For displaying purpose:
    pathstr_neuron = sprintf('%s\\tile%d_%d',newfolderPath,tileLabel,neuronLabel);
    display(pathstr_neuron)
end


end

