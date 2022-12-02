% PAHcalculateVesselArea.m
% Author: Erin McCaffrey
% Script reads in the vessel masks for all points, determines the area of
% all vessel regions and exports values as a csv for downstream stats.

%% Initiate paths to data and define points to analyze

path = ''; % general path
pathMask = ''; % location of vessel masks
resultsDir = path;
mkdir(resultsDir);
points =[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,23,24,25,...
    26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,...
    48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,...
    70,71,72,73,74,75,76,77,78,79,80,82,83,84,85,86,87,88,89,90];

%% Initiate an output matrix

vessel_area = []; 

%% Export and store vessel area stats 

for i=1:length(points)

    %load data
    disp(points(i));
    point = points(i);

    % load data for mask
    load([pathMask,'/Point',num2str(point),'/regionmasks.mat']);

    %get stats for objects in image
    stats = regionprops(filled_mask,'Area');

    %prepare output
    area_data = zeros(size(stats,1), 2);
    area_data(:,1) = point;
    area_data(:,2) = [stats(1:size(stats,1)).Area];
    
    %append
    vessel_area=[vessel_area;area_data];
end

%% Export as csv
colLabels = {'Point_num','area'};
TEXT.PnS = colLabels;
csvwrite_with_headers([resultsDir,'/vessel_area_data.csv'],vessel_area,TEXT.PnS)