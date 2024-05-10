% Neuron Harmonic Analysis given extracted neuron morphological file (.swc)
% and corresponding fluoresence images for individual neurons.
% Contact zhuang.xu@unsw.edu.au for further details

close all
set(0,'DefaultFigureVisible','on'); % turn 'on' or 'off' figure display

% coverslipID = "0759"; % coverslip number
coverslipID = "87"; % coverslip number
swcFolder = sprintf('D:\\Neuron_analysis\\neuronSWC%s',coverslipID); % double slash for sprintf function
markerFolder = sprintf('D:\\Neuron_analysis\\extractedNeurons%s',coverslipID); % Specify fluororescence image folder

swcFiles = dir(fullfile(swcFolder, '*.swc'));
swcName = string(natsortfiles({swcFiles.name}'));
swcpath = strcat(swcFolder, '\', swcName);

intfilesTau = dir(fullfile(markerFolder, '**\*Tau*.tiff'));
intfilesTubulin = dir(fullfile(markerFolder, '**\*tubOrig*.tiff')); % 16bit tubulin image
intfilesTubulinOutput = dir(fullfile(markerFolder, '**\*tubulin*.tiff')); % 8bit tubulin image
% intfilesTau = dir(fullfile(markerFolder, '**\*MAP2*.tiff'));

imgpathTau = fullfile(string(natsortfiles({intfilesTau.folder}')),'\',string(natsortfiles({intfilesTau.name}')));
imgpathTubulin= fullfile(string(natsortfiles({intfilesTubulin.folder}')),'\',string(natsortfiles({intfilesTubulin.name}')));
intfilesTubulinOutput = fullfile(string(natsortfiles({intfilesTubulinOutput.folder}')),'\',string(natsortfiles({intfilesTubulinOutput.name}')));% 8bit tubulin image

pixel2micron = 0.325; % micron/pixel size ratio for current microscope setup
somaRfacList = [1]; % soma radius scale factor for testing effects on spatial harmonics;
volfacList = [1]; % test scaling soma v;
areafacList = [1]; % test scaling soma surface area;

colNum = 30+6+8+12+1 ; % +6 for RTTList; +8 for 3 somaRfactor x 3 somavolumefactor, +1 otherAnkPeak
tableTau = {};
for iter = 1:length(swcFiles)
    fileID = iter;
    % if ~ismember(iter,keptFile)
    %     continue
    % end
    nodalDistList = [];
    close all
    
    % thicknessGC = 1;
    data = importswc(swcpath(fileID)); % data used for analysis
    data(:,3:4) = data(:,3:4) + 1; % all positions shifted by 1 as Matlab index starts from 1.
    insertNum = 5; % number of skeletal points inserted between two adjacent nodes.
    data = skelExpand(data,insertNum);
    dataIntExt = data; % data used for intensity extraction & AnkG peak detection;
    label = data(:,1)';
    type = data(:,2)';
    data_xyz = data(:,3:5);
    radius = data(:,6)';
    connct = data(:,7)';
    radiusp = radius';
    
    %% Extract intensity profile along the neuronal skeleton
    imgMarkerTau = imread(imgpathTau(fileID));
    imgMarkerTubulin = imread(imgpathTubulin(fileID));
    
    % Flip matrix vertically as the origin of the skeleton start from bottom left, while matrix origin starts from top left.
    imgMarkerTau = flip(imgMarkerTau,1);
    imgMarkerTubulin = flip(imgMarkerTubulin,1);
    
    imgMarkerTau = double(imgMarkerTau);
    imgMarkerTubulin = double(imgMarkerTubulin);
    
    x = round(data(:,3));
    y = round(data(:,4));
    
    figure(1)
    bw = zeros(size(imgMarkerTau));
    idxw = sub2ind(size(imgMarkerTau), y, x);
    bw(idxw) = 1;
    img8 = uint8(imgMarkerTubulin/max(imgMarkerTubulin(:))*256);
    img8_eq = adapthisteq(img8); % Enhance contrast for visualization using contrast-limited adaptive histogram equalization
    bw2 = imoverlay(img8_eq,bw,'r');
    imshow(flip(bw2))
    axis tight
    axis equal
    
    [intensityTau intensityTubulin] = extractSkelInt(dataIntExt,imgMarkerTau,imgMarkerTubulin);
    
    % [ydim xdim] = size(imgMarkerTau);
    x = pixel2micron*(1:size(imgMarkerTau,2));
    y = pixel2micron*(1:size(imgMarkerTau,1));
    [X, Y] = meshgrid(x,y);
    
    % For enhancing low intensity peaks
    % imgMarkerAnkG2 = imgMarkerAnkG;
    % imgMarkerAnkG2(imgMarkerAnkG2 > max(max(imgMarkerAnkG))/2) = max(max(imgMarkerAnkG))/1.5;
    Xplot = X;
    Yplot = Y;
    figure(2)
    saturationFac = 1;
    imgMarkerTauPlot = imgMarkerTau;
    imgMarkerTauPlot(imgMarkerTauPlot > saturationFac*max(imgMarkerTauPlot(:))) = saturationFac*max(imgMarkerTauPlot(:));
    surf(Xplot,Yplot,imgMarkerTauPlot,'edgecolor','none')
    axis tight
    axis equal
    colormap(jet)
    colorbar
    view(2)
    title('Tau Intensity','interpreter','latex')
    xlabel('X-coordinate /$\mu m$','interpreter','latex')
    ylabel('Y-coordinate /$\mu m$','interpreter','latex')
    
    l = sum(type == 1); %redundant in this case, keep as backup
    % dimension of neigbrs, assume each point has maximum of 3 neighbours (each process can only split into at most two branches)
    dim = length(label)-(l-1);
    D= ones(dim,1);
    dist = sparse(dim,dim);
    neibrs = zeros(dim,3);
    neibrs_soma = label(connct == 1);
    for i = label(l+1:end)
        nbr = [connct(i),label(connct == i)]; %hence neibrs(:,1) are parent cells
        neibrs(i,1:length(nbr)) = nbr;
    end
    
    for i = label(l+1:end)
        for j =  neibrs(i,neibrs(i,:)~=0)
            if j >= i
                continue
            elseif j == 1
                %             consider carefully here!
                %             dist(i,j) = norm(data_xyz(j,:) - data_xyz(i,:))- radius(1);
                dist(i,j) = norm(data_xyz(j,:) - data_xyz(i,:))/2;
            else
                dist(i,j) = norm(data_xyz(j,:) - data_xyz(i,:))/2;
            end
        end
    end
    dist = dist+dist';
    
    [containGC,idx_gc] = detectGC(fileID,data,dist,neibrs,neibrs_soma,markerFolder);
    
    %% Construct volumetric Laplacian and calculate eigenvectors
    for somaRfac = somaRfacList
        data = importswc(swcpath(fileID)); % data used for analysis
        data(:,3:4) = data(:,3:4) + 1; % all positions shifted by 1 as Matlab index starts from 1.
        data(1,6) = data(1,6)*somaRfac;
        insertNum = 5; %5, number of skeletal points inserted between two adjacent nodes.
        data = skelExpand(data,insertNum); % newly added       
        dataIntExt = data; % data used for intensity extraction & AnkG peak detection;
        label = data(:,1)';
        type = data(:,2)';
        data_xyz = data(:,3:5);
        radius = data(:,6)';
        connct = data(:,7)';
        l = sum(type == 1); %redundant in this case, keep as backup
        dim = length(label)-(l-1);
        D= ones(dim,1);
        for volfac = volfacList
            % preallocating volume and cross-sectional area at each skeletal point and related matrices
            volume = zeros(1,dim);
            Ai = pi*radius.^2;
            area = sparse(dim,dim);
            dist = sparse(dim,dim);
            lapl = sparse(dim,dim);
            neibrs = zeros(dim,3);%dimension of neigbrs = dim
            neibrs_soma = label(connct == 1);
            
            for i = label(l+1:end)
                nbr = [connct(i),label(connct == i)]; %hence neibrs(:,1) are parent cells
                neibrs(i,1:length(nbr)) = nbr;
            end
            
            % determine indices of all terminal points
            endpointList = [];
            for i = label(l+1:end)
                if sum(neibrs(i,:)~=0) == 1
                    endpointList = [endpointList i];
                else
                    continue
                end
            end
            
            for i = label(l+1:end)
                for j =  neibrs(i,neibrs(i,:)~=0)
                    if j >= i
                        continue
                    elseif j == 1
                        % different approximation for flux entering soma compartment
                        % dist(i,j) = norm(data_xyz(j,:) - data_xyz(i,:))- radius(1);
                        dist(i,j) = norm(data_xyz(j,:) - data_xyz(i,:))/2;
                        midRadius = radius(i);
                        area(i,j) = pi*midRadius^2;
                    else
                        dist(i,j) = norm(data_xyz(j,:) - data_xyz(i,:))/2;
                        midRadius = (radius(i)+radius(j))/2;
                        area(i,j) = pi*midRadius^2;
                    end
                end
            end
            dist = dist+dist';
            area = area+area';
            % calculate volume of soma, assuming a spherical geometry
            volume(1) = volfac*4/3*pi*radius(1)^3;
            
            for i = label(l+1:end)
                neibrsNot0 = neibrs(i,neibrs(i,:)~=0);%since 0 in the array means no neighbour
                if ismember(i,neibrs_soma) % check if it's connected to soma
                    neibrsNot01 = neibrsNot0(neibrsNot0~=1); %exclude the label for soma
                    if isempty(neibrsNot01) %if there is no neibour other than soma, return 1
                        %approximation of volume of points adjacent to soma
                        midpt = data_xyz(1,:) + radius(1)/(2*dist(i,1))*(data_xyz(i,:)-data_xyz(1,:));
                        volume(i-(l-1)) = 0.5*norm(midpt - data_xyz(i,:))*area(i,1);  % volume*0.5 to compensate the effect ghost point
                    else
                        midpt = data_xyz(1,:) + radius(1)/(2*dist(i,1))*(data_xyz(i,:)-data_xyz(1,:));
                        radius_midpoint = radius(i);
                        volume(i-(l-1)) = pi*norm(midpt - data_xyz(i,:))/3*(radius_midpoint^2+radius(i)^2+radius_midpoint*radius(i))+sum(1/3*dist(i,neibrsNot01).*...
                            (area(i,neibrsNot01)+Ai(i)+sqrt(Ai(i).*area(i,neibrsNot01))));
                    end
                    
                elseif length(neibrsNot0)==1 % if it's the end point
                    volume(i-(l-1)) = 0.5*1/3*dist(i,neibrsNot0).*(area(i,neibrsNot0)+Ai(i)...
                        +sqrt(Ai(i).*area(i,neibrsNot0))); % treat the point as the boundary, volume*0.5 to compensate the effect ghost point used to construct reflective boundary condition
                else
                    volume(i-(l-1)) = sum(1/3*dist(i,neibrsNot0).*(area(i,neibrsNot0)+Ai(i)...
                        +sqrt(Ai(i).*area(i,neibrsNot0))));
                end
            end
            
            % Construct adjacency matrix diagnol matrix and the volumetric Laplacian
            AM = (D./volume'.*(area.*spfun(@(x) 1./x, 2*dist)));
            DM = sparse(1:dim,1:dim, sum(AM,2),dim,dim);
            lapl = DM-AM;
            
            if length(label)<5000
                [eigvec eigval] = eig(full(lapl)); % calculate eigenpairs and sort eigenvectors based on magnitude of eigenvalues
                [d,ind] = sort(diag(eigval),'ascend');
                eigval_s = eigval(ind,ind);
                eigvec_s = eigvec(:,ind);
            else
                % use eigs to calculate eigenvector and eigenvalues if the matrix is too large
                %     opts.disp = 0;
                %     opts.isreal = 1;
                %     opts.maxit = 10000;
                %     opts.tol = 1e-14;
                %     opts.p = 500;
                [eigvec_s eigval_s] = eigs(lapl,12,'sm');
            end
            
            % change pixel to physical length units for analysis
            data(:,3:6) = pixel2micron*dataIntExt(:,3:6);
            data_xyz = data(:,3:5);
            radius = data(:,6)';
            somaR = data(data(:,7)==-1,6);
            
            % calculate path distance from every skeletal points to the estimated centre of soma
            pathDist = distTosoma(data_xyz, connct', type', label');
            zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); % determine the closest skeletal point to the nodal set
            idx_nodal = max(zci(eigvec_s(:,2))); % find index of the that point
            eigvec_s(:,1) =  sign(sum(eigvec_s(:,1)))*eigvec_s(:,1);
            eigvec_s(:,2) =  sign(eigvec_s(idx_nodal,2))*eigvec_s(:,2); % the point around nodal plane with larger label should have positive sign
            eigvalStore = diag(eigval_s)';
            
        end
        intensityTau_norm  = intensityTau/norm(intensityTau); % normalized Tau fluoresence for calculation of cross correlation
        linearDist = []; % Path dist by arraging individual arbors in an array
        for i = 1:length(label)
            if connct(i) == -1
                linearDist = [linearDist;0];
            else
                linearDist = [linearDist;linearDist(end)+norm(data(i,3:5)-data(connct(i),3:5))];
            end
        end
        
        weights = sum(full(dist),2); %
        idx_largeVariance = pathDist<1.2*radius(1);
        idx_fit = ~idx_largeVariance;
        radius_vol = radius';
        mseList = [];
        corrList = [];
        compList = [];
        paList = {};
        
        n =5;
        for k = 1:n
            X = radius_vol(idx_fit).*eigvec_s(idx_fit,1:k);
            [pa,stdx,mse] = lscov(X,intensityTau_norm(idx_fit),weights(idx_fit)); % pa -parameters for least-squre fitting
            compList(:,k) = X*pa;
            corrList(k) = weighted_corrcoef(intensityTau_norm(idx_fit),compList(:,k),weights(idx_fit));
            mseList(k) =  mse*(size(X,1)-size(X,2))/sum(weights);
            paList(k) = {pa};
        end
        corrList(n+1) = weighted_corrcoef(intensityTau_norm(idx_fit),compList(:,2)-compList(:,1),weights(idx_fit));
        
        pRrange = 2:length(intensityTau_norm);
        figure(3)
        plot(linearDist(pRrange),intensityTau_norm(pRrange),'b-','linewidth',0.8)
        hold on
        for k = 1:n
            plot(linearDist(idx_fit),compList(:,k),'linewidth',1)
        end
        hold off
        legend('data','0','1','2','3','4','5')
        
        corrList2 = [];
        comp0 = X(:,1)*pa(1);
        comp1 = X(:,1:2)*pa(1:2);
        comp2 = X(:,1:3)*pa(1:3);
        comp3 = X(:,1:4)*pa(1:4);
        comp4 = X(:,1:5)*pa(1:5);
        
        intensityTau_norm2 = intensityTau_norm(idx_fit) - X(:,1)*pa(1);
        corr0 = weighted_corrcoef(intensityTau_norm(idx_fit),comp0,weights(idx_fit));
        corr1 =  weighted_corrcoef(intensityTau_norm(idx_fit),comp1,weights(idx_fit)); % exclude region near soma for correlation calculation
        corr2 = weighted_corrcoef(intensityTau_norm(idx_fit),comp2,weights(idx_fit));
        corr3 = weighted_corrcoef(intensityTau_norm(idx_fit),comp3,weights(idx_fit));
        corr4 = weighted_corrcoef(intensityTau_norm(idx_fit),comp4,weights(idx_fit));
        corr0rm_1 = weighted_corrcoef(intensityTau_norm2,X(:,2),weights(idx_fit));
        corrList2 = [corr0 corr1 corr2 corr3 corr4 corr0rm_1];
        
        pRrange = 2:length(intensityTau_norm);
        
        figure(4)
        arborstart = find(data(:,7)==1);
        colors = {'#E13102','#0000CC','#00FF3F','#42F9F9','#E044A7','#FF00FF'};
        hold on
        idxBreakP = find((label-1) ~= connct);
        idxBreakP(1) = 2; % ignore soma
        comp0p = radius_vol.*eigvec_s(:,1)*pa(1);
        comp1p = radius_vol.*eigvec_s(:,1:2)*pa(1:2);
        for i = 1:length(idxBreakP)
            if i == length(idxBreakP)
                plot(linearDist(idxBreakP(i):end),intensityTau_norm(idxBreakP(i):end),'-','Color',colors{2},'linewidth',1.5);
                plot(linearDist(idxBreakP(i):end),comp0p(idxBreakP(i):end),'-','Color',colors{4},'linewidth',1.5)
                plot(linearDist(idxBreakP(i):end),comp1p(idxBreakP(i):end),'-','Color',colors{3},'linewidth',1.5)
            elseif i == 1 | (idxBreakP(i+1)-1) < arborstart(2)
                plot(linearDist([idxBreakP(i):idxBreakP(i+1)-1]),intensityTau_norm([idxBreakP(i):idxBreakP(i+1)-1]),'-','Color',colors{1},'linewidth',1.5);
                plot(linearDist([idxBreakP(i):idxBreakP(i+1)-1]),comp0p([idxBreakP(i):idxBreakP(i+1)-1]),'-','Color',colors{4},'linewidth',1.5)
                plot(linearDist([idxBreakP(i):idxBreakP(i+1)-1]),comp1p([idxBreakP(i):idxBreakP(i+1)-1]),'-','Color',colors{3},'linewidth',1.5)
            else
                plot(linearDist([idxBreakP(i):idxBreakP(i+1)-1]),intensityTau_norm([idxBreakP(i):idxBreakP(i+1)-1]),'-','Color',colors{2},'linewidth',1.5);
                plot(linearDist([idxBreakP(i):idxBreakP(i+1)-1]),comp0p([idxBreakP(i):idxBreakP(i+1)-1]),'-','Color',colors{4},'linewidth',1.5)
                plot(linearDist([idxBreakP(i):idxBreakP(i+1)-1]),comp1p([idxBreakP(i):idxBreakP(i+1)-1]),'-','Color',colors{3},'linewidth',1.5)
            end
        end
        hold off
        
        for k = 1:length(arborstart)
            xline(linearDist(arborstart(k)-1),'m--','LineWidth',1);
        end
        xlim([0 max(linearDist)+1])
        ylabel('Normalized intensity a.u.','interpreter','latex')
        xlabel('Stretched path distance ($\mu$m)','interpreter','latex')
        pdummy1 = line(nan,nan,'Linestyle','none','Marker','none','Color','none');
        pdummy2 = line(nan,nan,'Linestyle','none','Marker','none','Color','none');
        
        legend([pdummy1 pdummy2],['$r_0$ = ' num2str(round(corrList(1),3))],['$r_1$ = ' num2str(round(corrList(2),3))],...
            'Position',[0.1 0.82 0.2 0.05],'FontSize',12,'interpreter','latex')
        legend boxoff

        figure(5)
        idx_largeVariance = pathDist<1.2*radius(1); % only plot distal compartments
        % idx_largeVariance = pathDist<0; % show all compartments
        idx_distalProcesses = ~idx_largeVariance;
        saturationFac = 1;
        intensityTau_normPlot = intensityTau_norm;
        intensityTau_normPlot(intensityTau_normPlot > saturationFac*max(intensityTau_normPlot(idx_distalProcesses))) = saturationFac*max(intensityTau_normPlot(idx_distalProcesses));
        % % Plot Tau intensity profile
        scatter3(data(idx_distalProcesses,3),data(idx_distalProcesses,4),data(idx_distalProcesses,5),5,intensityTau_normPlot(idx_distalProcesses),'filled')
        
        colormap(jet)
        colorbar
        view(2)
        shading interp
        axis tight
        axis equal
        grid off
        xlabel('X ({\mum})')
        ylabel('Y ({\mum})')
        
        figure(6)
        % % Plot Tau intensity profile
        comp1Plot = comp1;
        saturationFac = 1;
        comp1Plot(comp1Plot > saturationFac*max(comp1Plot)) = saturationFac*max(comp1Plot);
        
        scatter3(data(idx_distalProcesses,3),data(idx_distalProcesses,4),data(idx_distalProcesses,5),5,comp1Plot,'filled')
        colormap(jet)
        colorbar
        view(2)
        shading interp
        axis tight
        axis equal
        grid off
        xlabel('X ({\mum})')
        ylabel('Y ({\mum})')
        
        figure(7)
        eigvecP = eigvec_s(:,2);
        n = 51; %give odd number for convenience
        [x,y,z] = sphere(n);
        somaColor = repmat(linspace(eigvecP(find(data(:,7)==-1)),min(eigvecP(find(data(:,7)==1))),(n+1)/2)',1,n+1);
        somaColorflp = flip(somaColor,1);
        somaColor = [somaColor; somaColorflp];
        subplot(1,4,1)
        scatter3(data(:,3),data(:,4),data(:,5),5,eigvecP,'filled')
        hold on
        % surf(somaR*x,somaR*y,somaR*z,ones(size(x))*eigvec_s(1,2))
        surf(somaR*x+data(1,3),somaR*y+data(1,4),somaR*z+data(1,5),somaColor)
        shading interp
        axis equal
        grid off
        colormap(jet)
        colorbar
        view([0 90])
        title('1st Inhomogeneous Harmonic','interpreter','latex')
        xlabel('X-coordinate /$\mu m$','interpreter','latex')
        ylabel('Y-coordinate /$\mu m$','interpreter','latex')
        
        %     display(datap(idx_nodal,:)) % display the row corresponding to the nodal plane
        
        subplot(1,4,2)
        eigvecP = eigvec_s(:,3);
        somaColor = repmat(linspace(eigvecP(find(data(:,7)==-1)),min(eigvecP(find(data(:,7)==1))),(n+1)/2)',1,n+1);
        somaColorflp = flip(somaColor,1);
        somaColor = [somaColor; somaColorflp];
        scatter3(data(:,3),data(:,4),data(:,5),5,eigvecP,'filled')
        hold on
        surf(somaR*x+data(1,3),somaR*y+data(1,4),somaR*z+data(1,5),somaColor)
        shading interp
        axis equal
        grid off
        colormap(jet)
        colorbar
        view(2)
        title('2nd Inhomogeneous Harmonic','interpreter','latex')
        xlabel('X-coordinate /$\mu m$','interpreter','latex')
        ylabel('Y-coordinate /$\mu m$','interpreter','latex')
        
        subplot(1,4,3)
        scatter3(data(:,3),data(:,4),data(:,5),5,radiusp.*eigvec_s(:,2),'filled')
        colormap(jet)
        colorbar
        view(2)
        shading interp
        axis tight
        axis equal
        grid off
        title('1st Harmonic Intensity')
        
        idx_largeVariance = pathDist<1*radius(1);
        idx_distalProcesses = ~idx_largeVariance;
        % Plot Tau intensity profile
        subplot(1,4,4)
        scatter3(data(idx_distalProcesses,3),data(idx_distalProcesses,4),data(idx_distalProcesses,5),5,intensityTau(idx_distalProcesses),'filled')
        colormap(jet)
        colorbar
        view(2)
        shading interp
        axis tight
        axis equal
        grid off
        xlabel('X /{\mum}')
        ylabel('Y /{\mum}')
        title('Tau Intensity')
          
        figure(8)
        img8 = uint8(imgMarkerTubulin/max(imgMarkerTubulin(:))*256);
        img8_eq = adapthisteq(img8); % Enhance contrast for visualization using contrast-limited adaptive histogram equalization
        output3 = imshow(flip(img8_eq));
        axis tight
        axis equal 
        
        % determine maximum length of each arbor  
        arborNum = length(neibrs_soma);
        arborList = []; % The indices of each distinct branch is stored in rows of arborList
        for i = 1:arborNum
            if i == arborNum
                arborList(i,1:(dim-neibrs_soma(i)+1)) = neibrs_soma(i):dim; % This makes use of of the structure of swc files
            else
                arborList(i,1:(neibrs_soma(i+1)-neibrs_soma(i))) = neibrs_soma(i):(neibrs_soma(i+1)-1);
            end
        end
        
        % Find the longest path length (presumably dendritic) that doesn't go through nodal plane
        pathLengthArbors = [];
        idx_longArborEnd = NaN;
        arborLabels = [];
        for j = 1:arborNum
            arborLabels = arborList(j,arborList(j,:)~=0);
            arborLabels = arborLabels(arborLabels~=0);
            pathLengthArbors(j) = max(pathDist(arborLabels));
        end
        arborLsort =  sort(pathLengthArbors,'descend');
        neuronL = sum(arborLsort(1:2));
        figure(9)
        title(['$L_{neuron}$ = ' num2str(round(neuronL)) '$\mu$m'],'interpreter','latex')

        % Output figures, save in fig and svg format
        outFolder = 'C:\Users\ericx\Documents\NeuronHarmonicAnalysis\Examples'; % specify folder for figure output
        coverslipNeuron = 'exampleNeuron1'; % name for the neuron analyzed
        
        figure(10)
        fileName = [coverslipNeuron '1D Tau fitting'];
        filePath = fullfile(outFolder,fileName);
        saveas(gcf,filePath,'fig')
        saveas(gcf,filePath,'svg')
        
        figure(11)
        fileName = [coverslipNeuron 'Tau image'];
        filePath = fullfile(outFolder,fileName);
        saveas(gcf,filePath,'fig')
        saveas(gcf,filePath,'svg')  
        
        figure(12)
        fileName = [coverslipNeuron 'Tau skeletal intensity'];
        filePath = fullfile(outFolder,fileName);
        saveas(gcf,filePath,'fig')
        saveas(gcf,filePath,'svg')
        
        figure(13)
        fileName = [coverslipNeuron 'harmonics skeletal intensity'];
        filePath = fullfile(outFolder,fileName);
        saveas(gcf,filePath,'fig')
        saveas(gcf,filePath,'svg')
        
        figure(14)
        fileName = [coverslipNeuron 'Tubulin image_2'];
        filePath = fullfile(outFolder,fileName);
        saveas(gcf,filePath,'fig')
        saveas(gcf,filePath,'svg')
        
        fileName = [coverslipNeuron 'Tubulin image.tiff'];
        outputFullFileName = fullfile(outFolder,fileName);
        inputFullFileName = intfilesTubulinOutput(iter);
        copyfile(inputFullFileName, outputFullFileName);
        % store table of analyzed data
        tableTau(iter,:) =  {eigvalStore(1:5),eigvec_s(:,1:5),radius',weights,pathDist,pathLengthArbors,compList,paList,mseList,corrList,corrList2,containGC};
        fprintf('iter: %d , fileID: %d\n',iter,fileID)
    end
end
