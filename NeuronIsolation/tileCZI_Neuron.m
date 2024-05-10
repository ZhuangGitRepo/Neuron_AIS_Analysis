% This code tiles the original CZI images

coverslipNumber = '87';
path_to_czi = 'C:\NeuronCzi\TF_coverslips\2019_11_22_87.czi';
imageName = '2019_11_22_87';

% create the folder for tiled neuron images and extracted individual
% neurons
d = sprintf('D:\\\\Neuron_analysis\\\\tiledNeurons%s',coverslipNumber);
mkdir(d)
d2 = sprintf('D:\\\\Neuron_analysis\\\\extractedNeurons%s',coverslipNumber);
mkdir(d2)
% To extract neuron using vaa3D and output extracted neuron via cmd, use single backslash in
% forder path
d3 = sprintf('D:\\Neuron_analysis\\neuronSWC%s',coverslipNumber);
mkdir(d3)
d4 = sprintf('D:\\Neuron_analysis\\extractedNeurons%s',coverslipNumber); % Format for cmd input

% Below shows the fluorophore marker order.
% 2019_11_22_87 - 3 DIV - Staining combination 1 [¦Â3-tubulin-A488; AnkG-A647; Tau1-A555]

% Marker sequance for the coverslip, ordered in ascending wavelengths
MarkerStrings = ["DAPI","tubulin","Tau","AnkG"];


for k = 1:3% set timeSeries Number

    autostich = true;
    attachments = true;
    flattenedResolutions = false;
    [r seriesCount] = openCZI(path_to_czi,autostich,flattenedResolutions,attachments);  % Seriescount: in this case the number of different types of fluorophhores 

    timeSeries = k;  % select 1st, 2nd or 3rd 
    r.setSeries(timeSeries-1)
    resolutionCount = r.getResolutionCount();
    planeCount = r.getImageCount(); % Number of image in the ciz file
    fileID = 1;
    fprintf(fileID, 'Series: %g\n', timeSeries);
    fprintf(fileID, 'Number of fluorescence microscope images: %g\n', planeCount);
    fprintf(fileID, 'Number of resolution levels: %g\n', resolutionCount);
    fprintf(fileID, 'stackSizeX: %g\n', r.getSizeX());
    fprintf(fileID, 'stackSizeY: %g\n', r.getSizeY());
    fprintf(fileID, 'stackSizeZ: %g\n', r.getSizeZ());
    fprintf(fileID, 'stackSizeC: %g\n', r.getSizeC());
    fprintf(fileID, 'stackSizeT: %g\n', r.getSizeT());

    %access physical voxel and stack sizes of the data
    omeMeta = r.getMetadataStore();
    stackSizeX = omeMeta.getPixelsSizeX(timeSeries-1).getValue(); % image width, pixels
    stackSizeY = omeMeta.getPixelsSizeY(timeSeries-1).getValue(); % image height, pixels
    stackSizeZ = omeMeta.getPixelsSizeZ(timeSeries-1).getValue(); % number of Z slices
    stackSizeC = omeMeta.getPixelsSizeC(timeSeries-1).getValue(); % number of colors (fluorescence markers)
    stackSizeT = omeMeta.getPixelsSizeZ(timeSeries-1).getValue(); % number of time sequence in the plane
    voxelSizeXdefaultValue = omeMeta.getPixelsPhysicalSizeX(timeSeries-1).value();           % returns value in default unit
    voxelSizeXdefaultUnit = omeMeta.getPixelsPhysicalSizeX(timeSeries-1).unit().getSymbol(); % returns the default unit type
    voxelSizeX = omeMeta.getPixelsPhysicalSizeX(timeSeries-1).value(ome.units.UNITS.MICROMETER); % in ?m
    voxelSizeXdouble = voxelSizeX.doubleValue();                                  % The numeric value represented by this object after conversion to type double
    voxelSizeY = omeMeta.getPixelsPhysicalSizeY(timeSeries-1).value(ome.units.UNITS.MICROMETER); % in ?m
    voxelSizeYdouble = voxelSizeY.doubleValue();                                  % The numeric value represented by this object after conversion to type double
    if omeMeta.getPixelsPhysicalSizeZ(timeSeries-1)~=[]
        voxelSizeZ = omeMeta.getPixelsPhysicalSizeZ(timeSeries-1).value(ome.units.UNITS.MICROMETER); % in ?m
        voxelSizeZdouble = voxelSizeZ.doubleValue();                                  % The numeric value represented by this object after conversion to type double
    else
        voxelSizeZdouble = 0; %No z dimension
    end

    resolutionLevel = 0; %resolution level ranges from 0 to (resolutionCount-1), choose 0 for highest resolution
    r.setResolution(resolutionLevel);

    tilesPerRow = 5;
    tilesPerCol = 5;

    width = floor(stackSizeX/tilesPerRow);
    height = floor(stackSizeY/tilesPerCol);

    tileNum = 1;

    for j = 1:height:stackSizeY+1 
        for i = 1:width: stackSizeX+1 
            % Tile 2019_11_22_86
            if i+width > (stackSizeX+1) || j+height > (stackSizeY+1) 
                continue
            end
            % Marker sequance for 2019_11_22_87.czi
            imgSeries_1 = bfGetPlane(r,1, i,j,width,height); %1 for DAPI markers
            imgSeries_2 = bfGetPlane(r,2, i,j,width,height); %2 for tubulin markers
            imgSeries_3 = bfGetPlane(r,3, i,j,width,height); %3 for Tau markers
            imgSeries_4 = bfGetPlane(r,4, i,j,width,height); %4 for AnkG markers
            tiledNeuron = cat(3,imgSeries_1,imgSeries_2,imgSeries_3,imgSeries_4);

            tileName = sprintf('%s_time%d_tile%d.tiff',imageName,timeSeries,tileNum);
            imwrite(tiledNeuron,fullfile(d,tileName));
            tileNum = tileNum + 1;
        end
    end
end


%% Extract neurons from tiled images
voxelSizeYdouble = 0.325;
voxelSizeXdouble = 0.325;

tiles = dir(fullfile(d, '*.tiff'));
NumExtracted = [];
radiusPixelList = [];
areaPixelList = []; 
DAPIPixelList = [];
circularityList = [];
blurrinessList = [];
contactingNeuron_total = 0;
unidentifiedObj_total = 0;
loopNum_total = 0;
roundCell_total = 0;
darkRegion_total = 0;
blurImage_total = 0;

neuronLabel = 0; % Label neuron in each tiled image
tileNames = {tiles.name};
tileFolder = {tiles.folder};
tiles = strcat(tileFolder, '\', tileNames);
tiles = natsortfiles(tiles);

for i = 1:length(tiles)
    tileLabel = i;
    % create neuron storage folder for each tile
    newfolder = sprintf('tile%d',i);
    newfolderPath = fullfile(d2,newfolder);
    mkdir(newfolderPath);
    region2D = [];
    imagePath = tiles{i};
    multiMarkersIMG = imread(imagePath);
    
    % Specify intensity threshold for pixels identified as part of a neuron, value ranges from 0-255,
    %high value would indicate noisy background signals 
%     intThresh = 1; %8bit:1-5
    intThresh = 1300;
    % Specify smallest objects classified as neurons, given in pixels 
    pixelThreshL = 10000;
    % ExtendLenth - extend the window for each neuron to make sure all signal is
    % covered, given in pixels
    extendLenth = 20;
    
    % detect soma by its physical size, specify the smallest neucleus size by nucleusSizeThreshL
    % General neucleus size of 3DIV neuron ranges between 700-1100 pixels
%     nucleusSizeThreshL = 7; % minimum nucleus size in mircon
    nucleusSizeThreshL = 5; % minimum nucleus size in mircon
    nucleusPixelThreshL =  pi*(nucleusSizeThreshL/voxelSizeYdouble/2)*(nucleusSizeThreshL/voxelSizeXdouble/2);
    nucleusPixelThreshL = round(nucleusPixelThreshL);
    
%     nucleusSizeThreshH = 14; % maximum nucleus size 
    nucleusSizeThreshH = 20; % maximum nucleus size 
    nucleusPixelThreshH = pi*(nucleusSizeThreshH/voxelSizeYdouble/2)*(nucleusSizeThreshH/voxelSizeXdouble/2); % Area of ellipse pi*a*b
    nucleusPixelThreshH = round(nucleusPixelThreshH);
  

    [radiusPixel areaPixel DAPIPixel region2D circularity blurriness contactingNeuron unidentifiedObj roundCell loopNum blurImage darkRegion] = detectNeuron(multiMarkersIMG,intThresh,pixelThreshL,...
    pixelThreshH,extendLenth,nucleusPixelThreshL,nucleusPixelThreshH,newfolderPath,neuronLabel,tileLabel,MarkerStrings);

    radiusPixelList = [radiusPixelList; radiusPixel];
    areaPixelList = [areaPixelList; areaPixel]; 
    DAPIPixelList = [DAPIPixelList; DAPIPixel];
    NumExtracted = [NumExtracted; size(region2D,1)];
    circularityList = [circularityList; circularity];
    blurrinessList = [blurrinessList; blurriness];
    
    contactingNeuron_total = contactingNeuron_total + contactingNeuron;
    unidentifiedObj_total = unidentifiedObj_total + unidentifiedObj;
    roundCell_total = roundCell_total + roundCell;
    loopNum_total = loopNum_total + loopNum;
    darkRegion_total = darkRegion_total + darkRegion;
    blurImage_total = blurImage_total + blurImage;
end

dataLists = [radiusPixelList areaPixelList DAPIPixelList circularityList blurrinessList];
dataNum = [contactingNeuron_total unidentifiedObj_total roundCell_total loopNum_total darkRegion_total blurImage_total];


%% Automated tracing of extracted neurons (APP2)

files_neuron = dir(fullfile(d4,'**\*tubulin*.tiff'));
neuron_name = natsortfiles({files_neuron.name}');
neuron_folder = natsortfiles({files_neuron.folder}');

for i = 1:length(files_neuron)
    imagepath = string(fullfile(neuron_folder(i), neuron_name(i)));
    suffix = '.tiff';
    swcname = string(strrep(neuron_name(i),suffix,''));
    swcname = sprintf('%s_neuron%d.swc',swcname,i);
    vaa3dpath = 'C:\Vaa3D_V3.601_Windows_MSVC_64bit';
    functionpath = 'C:\Vaa3D_V3.601_Windows_MSVC_64bit\plugins\neuron_tracing\Vaa3D_Neuron2\vn2.dll';
    inpath = imagepath;
    outpath = sprintf('%s\\%s',d3,swcname);
    %vaa3d app2 function input format: vaa3d -x plugin_name -f app2 -i <inimg_file> -o <outswc_file>...
    % -p <inmarker_file> <channel> <bkg_thresh> <b_256cube> <b_RadiusFrom2D> ...
    % <is_gsdt> <is_gap> <length_thresh> <is_resample> <is_brightfield> <is_high_intensity>
    % do not downsample b_256cube(downsample) but resample, radius from 2D,with gsdt, no gap, length_threshold 0.5.
    
    % In the cmd, cd to vaa3dpath, then type vaa3d_msvc /x vn2 /f help to
    % see info about APP2/APP1
    parameters = sprintf('NULL 0 AUTO 0 1 1 0 3 1 0 0');
    command = sprintf('cd %s && vaa3d_msvc /x %s /f app2 /i %s /o %s /p %s',vaa3dpath,functionpath,inpath,outpath,parameters);
    [status cmdout] = system(command);
end



