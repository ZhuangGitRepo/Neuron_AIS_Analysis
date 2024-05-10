function [containGC,idx_gc] = detectGC(fileID,data,dist,neibrs,neibrs_soma,d)
% This function detect if the neuron contains noticeable growth cone region from the SWC file
% by setting a threshold for density of skeletal points in  the region.
% neuron contains a growth cone need to be checked manually for quality control
idx_gc = [];
containGC = 0;
label = data(:,1)';
type = data(:,2)';
data_xyz = data(:,3:5);
connct = data(:,7)';
l = sum(type==1);
dim = length(label)-(l-1);

% Import dilated tubulin mask
intfiles = dir(fullfile(d, '**\*overlap*.tiff'));
intname = string(natsortfiles({intfiles.name}'));
intfolder = string(natsortfiles({intfiles.folder}'));
imgpath = fullfile(intfolder(fileID),'\',intname(fileID));

neuronMask = imread(imgpath);
[ydim xdim] = size(neuronMask);
neuronMask = logical(neuronMask(:,:,2)+neuronMask(:,:,3));
figure(100)
subplot(1,2,1)
imshow(neuronMask)
r = 7;
SE = strel('disk',r);
neuronMaskErode = imerode(neuronMask,SE);
D = bwdist(~neuronMaskErode);
RGB = repmat(rescale(D), [1 1 3]);
subplot(1,2,2)
imshow(RGB)
hold on, imcontour(D)
[C h] = contour(D);
D1 = D;
% contourThresh = median(h.LevelList)*7/5; % Original setting
contourThresh = median(h.LevelList)*5/5; % changed 2022/05/11

D1(D1<contourThresh)=false;

r = floor(contourThresh);
SE = strel('disk',r);
D2 = imdilate(D1,SE);
r = 14;
SE = strel('disk',r);
D3 = imdilate(D2,SE);

[B,L] = bwboundaries(D3,'noholes');

figure(101)
imshow(label2rgb(L, @jet, [.5 .5 .5]))

for k = 1:length(B)
    fatRegion = B{k};
    %    plot(fatRegion(:,2), ydim-fatRegion(:,1), 'r', 'LineWidth', 2)
end
% plot(data(:,3),data(:,4),'k.')
% axis equal

xq = data(:,3);
yq = data(:,4);
for i = 1:length(B)
    fatRegion = B{i};
    rng default
    [in,on] = inpolygon(xq,yq,fatRegion(:,2),ydim-fatRegion(:,1));
    fatArea = polyarea(fatRegion(:,2),ydim-fatRegion(:,1));
    containSoma = inpolygon(data(1,3),data(1,4),fatRegion(:,2),ydim-fatRegion(:,1));
    if containSoma == true
        containGC = 0;
        idx_gc = [];
        %         display('contain soma')
        continue
    end
    pathDist = distTosoma(data_xyz, connct', type', label');
    
    endpointList = [];
    for i = label(l+1:end)
        if sum(neibrs(i,:)~=0) == 1
            endpointList = [endpointList i];
        else
            continue
        end
    end
    arborNum = length(neibrs_soma);
    arborList = [];
    for i = 1:arborNum
        if i == arborNum
            arborList(i,1:(dim-neibrs_soma(i)+1)) = neibrs_soma(i):dim;
        else
            arborList(i,1:(neibrs_soma(i+1)-neibrs_soma(i))) = neibrs_soma(i):(neibrs_soma(i+1)-1);
        end
    end
    idx_in = find(in==1);
    
    % determine which arbor is the growth-cone-like region located
    onArborList = [];
    onArborPortion = []; %percentage of the points in the growth-cone-like region on certain arbor
    onArbor = [];
    for i = 1:arborNum
        if sum(ismember(idx_in,arborList(i,:))) ~=0
            onArborList = [onArborList i];
            onArborPortion = [onArborPortion sum(ismember(idx_in,arborList(i,:)))];
        end
    end
    
    onArbor = onArborList(onArborPortion==max(onArborPortion)); % Find the main arbor where the growth-cone-like region is located
    endpoint_pathDist = [];
    endpointViaFlatRegion = []; % find all end point that go throught the flatRegion
    for i = endpointList
        itemp = i;
        if ismember(i,arborList(onArbor,:)) == 0
            continue
        elseif ismember(i,idx_in) == 1
            endpointViaFlatRegion = [endpointViaFlatRegion i];
            continue
        else
            runIter = 1;
            while(runIter)
                itemp = connct(itemp);
                if itemp == -1
                    runIter = 0;
                elseif ismember(itemp,idx_in) == 1
                    endpointViaFlatRegion = [endpointViaFlatRegion i];
                    runIter = 0;
                else
                    continue
                end
            end
        end
    end
    
    if isempty(endpointViaFlatRegion)
        display("Morphology is not properly constructed")
        continue
    end
    endpoint_Longest = endpointViaFlatRegion(pathDist(endpointViaFlatRegion)==max(pathDist(endpointViaFlatRegion)));
    longestPathDistViaFlatRegion = pathDist(endpoint_Longest);
    
    if min(pathDist(in))/longestPathDistViaFlatRegion > 0.4 && max(pathDist(in))/longestPathDistViaFlatRegion > 0.85 % changed at 2022/05/11
        idx_gc = inpolygon(xq,yq,fatRegion(:,2),ydim-fatRegion(:,1)); % Indices of growth cone to be modified
        containGC = 1;
        break
    elseif length(B) == 1
        containGC = 0;
        idx_gc = [];
    end
end

end

