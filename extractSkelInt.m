function [skelInt1 skelIntTubulin] = extractSkelInt(data,imgMarker1,imgMarkerTubulin)
label = data(:,1)';
type = data(:,2)';
data_xyz = data(:,3:5);
radius = data(:,6)';
connct = data(:,7)';
dim = length(label);
D= ones(dim,1);
normVecList = zeros(dim,3,3); % normalized direction vector list. Dimension1: label; 2: neighbours; 3:xyz

% preallocating volume and cross section area at each point and related matrices
volume = zeros(1,dim);
areaSide = zeros(1,dim);

Ai = pi*radius.^2;
area = sparse(dim,dim);
circum = sparse(dim,dim);
dist = sparse(dim,dim);
distSurface =  sparse(dim,dim);
neibrs_soma = label(connct == 1);
neibrs = zeros(dim,3);% dimension of neigbrs = dim, assume each point has maximum of 3 neighbours (each process can only split into at most two branches)

for i = label(2:end)
    nbr = [connct(i),label(connct == i)]; %hence neibrs(:,1) are parent cells
    neibrs(i,1:length(nbr)) = nbr;
end

% Determine index of all terminal points
endpointList = [];
for i = label(2:end)
    if sum(neibrs(i,:)~=0) == 1
        endpointList = [endpointList i];
    else
        continue
    end
end

for i = label(2:end)
    for j =  neibrs(i,neibrs(i,:)~=0)
        if j >= i
            continue
        elseif j == 1
            dist(i,j) = norm(data_xyz(j,:) - data_xyz(i,:))/2;
            midRadius = radius(i);
            area(i,j) = pi*midRadius^2;
            circum(i,j) = 2*pi*midRadius;
            normVecList(i,j,:) = (data_xyz(j,:) - data_xyz(i,:))/norm(data_xyz(j,:) - data_xyz(i,:));
            normVecList(j,i,:) = -normVecList(i,j,:);
        else
            dist(i,j) = norm(data_xyz(j,:) - data_xyz(i,:))/2;
            midRadius = (radius(i)+radius(j))/2;
            area(i,j) = pi*midRadius^2;
            circum(i,j) = 2*pi*midRadius;
            normVecList(i,j,:) = (data_xyz(j,:) - data_xyz(i,:))/norm(data_xyz(j,:) - data_xyz(i,:));
            normVecList(j,i,:) = -normVecList(i,j,:);
        end
    end
end
dist = dist+dist';
area = area+area';
circum=circum+circum';
midrad = circum/(2*pi);

areaSide(1) = 4*pi*radius(connct==-1)^2;
volume(1) = 4/3*pi*radius(1)^3;

for i = label(2:end)
    neibrsNot0 = neibrs(i,neibrs(i,:)~=0);%since 0 in the array means no neighbour
    if ismember(i,neibrs_soma) %check if it's connected to soma
        neibrsNot01 = neibrsNot0(neibrsNot0~=1); %exclude the label for soma
        if isempty(neibrsNot01) %if there is no neibour other than soma, return 1
            % approximation of volumes of points adjacent to soma
            midpt = data_xyz(1,:) + radius(1)/(2*dist(i,1))*(data_xyz(i,:)-data_xyz(1,:));
            volume(i) = norm(midpt - data_xyz(i,:))*area(i,1);
            areaSide(i) = norm(midpt - data_xyz(i,:))*circum(i,1);
        else
            midpt = data_xyz(1,:) + radius(1)/(2*dist(i,1))*(data_xyz(i,:)-data_xyz(1,:));
            radius_midpoint = radius(i);
            volume(i) = pi*norm(midpt - data_xyz(i,:))/3*(radius_midpoint^2+radius(i)^2+radius_midpoint*radius(i))+sum(1/3*dist(i,neibrsNot01).*...
                (area(i,neibrsNot01)+Ai(i)+sqrt(Ai(i).*area(i,neibrsNot01))));
            areaSide(i) = norm(midpt - data_xyz(i,:))*circum(i,1) +sum(pi*(radius(i)+...
                midrad(i,neibrsNot01)).*sqrt((radius(i)-midrad(i,neibrsNot01)).^2+(dist(i,neibrsNot01)).^2));
        end
    elseif length(neibrsNot0)==1 % If the point is the end point
        volume(i) = 1/3*dist(i,neibrsNot0).*(area(i,neibrsNot0)+Ai(i)...
            +sqrt(Ai(i).*area(i,neibrsNot0)));
        areaSide(i) = pi*(radius(i)+midrad(i,neibrsNot0)).*sqrt((radius(i)-...
            midrad(i,neibrsNot0)).^2+(dist(i,neibrsNot0)).^2);
    else  % all other non-branching, non-end points
        volume(i) = sum(1/3*dist(i,neibrsNot0).*(area(i,neibrsNot0)+Ai(i)...
            +sqrt(Ai(i).*area(i,neibrsNot0))));
        areaSide(i) = sum(pi*(radius(i)+midrad(i,neibrsNot0)).*sqrt((radius(i)...
            -midrad(i,neibrsNot0)).^2+(dist(i,neibrsNot0)).^2));
    end
end

skelInt1 = zeros(dim,1);
skelIntTubulin = zeros(dim,1);

x = round(data(:,3));
y = round(data(:,4));
BW_skeleton = zeros(size(imgMarker1));
idxw = sub2ind(size(imgMarker1), y, x);
BW_skeleton(idxw) = 1;
img8 = uint8(imgMarkerTubulin/max(imgMarkerTubulin(:))*256);
img8_eq = adapthisteq(img8); % Enhance contrast for visualization using contrast-limited adaptive histogram equalization
mask_overlay = imoverlay(img8_eq,BW_skeleton,'r');
% figure()
% imshow(flip(mask_overlay))
% axis tight
% axis equal

radiusExt = radius;
radiusExt(2:end) = round(radius(2:end)*1.5); % width of extracted intensity region

for i = label(1:end)
    % Only find directional vector pointing towards child points except the
    % end points
    if i == 1
        childNeibrs = neibrs_soma;
    elseif length(neibrs(i,neibrs(i,:)~=0)) == 1
        childNeibrs = neibrs(i,neibrs(i,:)~=0);
    else
        childNeibrs = neibrs(i,neibrs(i,:)~=0 & neibrs(i,:)>i);
    end
    
    for j = childNeibrs
        v = squeeze(normVecList(i,j,1:2))';
        vperp = [-v(2), v(1)];

        coords_i = data_xyz(i,1:2);
        coords_perpSegEnds = [coords_i-vperp*radiusExt(i); coords_i+vperp*radiusExt(i)];
        x_perpSegEnds = round(coords_perpSegEnds(:,1));
        y_perpSegEnds = round(coords_perpSegEnds(:,2));
        [x_perpSeg y_perpSeg] = bresenham(x_perpSegEnds(1),y_perpSegEnds(1),x_perpSegEnds(2),y_perpSegEnds(2));
        mask = false(size(imgMarker1));
        mask(sub2ind(size(mask), y_perpSeg, x_perpSeg)) = 1;
        %         mask_overlay = imoverlay(mask_overlay,mask); % For visualization of the compartmentalization
        if i == 1 | length(x_perpSegEnds)<5;
            coords_i = data_xyz(i,1:2);
            coords_perpSegEnds = [coords_i-vperp*radius(i); coords_i+vperp*radius(i)];
            x_perpSegEnds = round(coords_perpSegEnds(:,1));
            y_perpSegEnds = round(coords_perpSegEnds(:,2));
            [x_perpSeg y_perpSeg] = bresenham(x_perpSegEnds(1),y_perpSegEnds(1),x_perpSegEnds(2),y_perpSegEnds(2));
            mask = false(size(imgMarker1));
            mask(sub2ind(size(mask), y_perpSeg, x_perpSeg)) = 1;
            skelInt1(i) = skelInt1(i) + max(imgMarker1(mask));
            skelIntTubulin(i) = skelIntTubulin(i) + max(imgMarkerTubulin(mask));
        else
            meanInt = estimateRange(imgMarker1,imgMarkerTubulin,i);
            skelInt1(i) = skelInt1(i) + meanInt(1);
            skelIntTubulin(i) = skelIntTubulin(i) + meanInt(3);
        end
        mask_overlay = imoverlay(mask_overlay,mask);
        %         figure(10000)
        %         imshow(mask_overlay)
    end
    if length(childNeibrs) > 1
        skelInt1(i) = skelInt1(i)/length(childNeibrs);
        skelIntTubulin(i) = skelIntTubulin(i)/length(childNeibrs);
    end
end

function meanInt = estimateRange(imgMarker1,imgMarkerTubulin,i)
temp = sub2ind(size(mask), y_perpSeg, x_perpSeg);
    for k = 1:2 % three marker files
        if k == 1
            tempInt = imgMarker1(temp)';
        else
            tempInt = imgMarkerTubulin(temp)';
        end
        tempDist = [0 sqrt((x_perpSeg(2:end)-x_perpSeg(1)).^2+(y_perpSeg(2:end)-y_perpSeg(1)).^2)'];
        height0 = max(tempInt) - min(tempInt);
        peakloc0 = tempDist(1)+radiusExt(i);
        std0 = radius(i);
        bkg0 = min(tempInt);
        gaussfitfun = @(c,tempDist) c(1)*exp(-((tempDist-c(2))/(sqrt(2)*c(3))).^2) + c(4);
        lb = [0 0 0 0];
        ub = [1.2*max(tempInt) tempDist(end) tempDist(end) bkg0];
        c0 = [height0 peakloc0 std0 bkg0];
        opts = optimset('Display','off');
        c = lsqcurvefit(gaussfitfun,c0,tempDist,tempInt,lb,ub,opts);
        distq = linspace(tempDist(1)-radiusExt(i),tempDist(end)+radiusExt(i),200);
        fittedGaussian =  gaussfitfun(c,distq);
        FWHM = 2*sqrt(2*log(3/2))*c(3);  % Full width at 2/3 maximum for gaussian fit
        extRange = distq>=(c(2)-FWHM/2) & distq<=(c(2)+FWHM/2);
        meanInt(k) = max(fittedGaussian(extRange));
    end
end

end




