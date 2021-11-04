% MIBItileOneMarkerAcrossPoints
% Script creates a tiles images of each channel in different points in the
% dimensions specified by the user. All points are scaled the same.

corePath = {'SampleData/extracted/Point1/dataNoBg.mat', ...
    'SampleData/extracted/Point2/dataNoBg.mat'}; % cores to work on. Can add several paths, separated by commas.
massFile = 'SampleData/SamplePanel.csv'; % panel csv
xTileNum = 1; % Number of rows in tile
yTileNum = 2; % Number of columns in tile
outDir = 'SampleData/extracted/TiledImages';
defaultCap = 5; % Cap to use if no other cap is specified in the massFile
xSize = 1030; % X-Size of the largest image to be tiled. Can add a few pixels to generate a border
ySize = 1030; % Y-Size of the largest image to be tiled. Can add a few pixels to generate a border

%% script

massDS = MibiReadMassData(massFile);
coreNum = length(corePath);
mkdir (outDir);
% load all cores
p=cell(coreNum,1);
for i=1:coreNum
    p{i} = load([corePath{i}]);
end

for i=1:length(massDS)
    % check if massDS has the cap variable and if it isn't empty
    if ismember('Cap', massDS.Properties.VarNames) && ~isempty(massDS.Cap(i))
        currCap = massDS.Cap(i);
    else
        currCap = defaultCap;
    end
    % Generate the data to plot
    tiledIm = zeros(xTileNum*xSize,yTileNum*ySize);
    for j=1:coreNum
        % cap and pad dat
        data = p{j}.countsNoBg(:,:,i);
        data(data>currCap)=currCap;
        % if data is smaller than expected, pad it
        dataPad = zeros(xSize,ySize);
        dataPad([1:size(data,1)],[1:size(data,2)]) = data;
        % get position
        xpos = floor((j-1)/yTileNum)+1;
        ypos = mod(j,yTileNum);
        if (ypos == 0)
            ypos = yTileNum;
        end
        tiledIm([(xpos-1)*xSize+1:(xpos-1)*xSize+xSize],[(ypos-1)*ySize+1:(ypos-1)*ySize+ySize])=dataPad;
    end
    imwrite(uint16(tiledIm),[outDir,'/',massDS.Label{i},'_tiled.tif']);
end
