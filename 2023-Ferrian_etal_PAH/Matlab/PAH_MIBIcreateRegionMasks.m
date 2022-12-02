% MIBIcreateGranulomaRegionMasks.m
% Date created: 190620
% Author: Erin McCaffrey 

% Reads in a mask tiff generated externally, allows modification to smooth 
% mask, creates a peripheral zone mask, and then exports as tiff and
% matlab variables

%define paths and points to analyse
%path = '/Volumes/GoogleDrive/My Drive/Angelo Lab/MATLAB/MIBI_SpatialAnalysis/masks/';
path = '/Users/sferrian/Desktop/FinalData/PAHanalyses/AllPoints/';
points = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,82,83,84,85,86,87,88,89,90]; 
%points = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]; %cohort data
imSize=1024;

%define masking parameters (these are tunable)
se1 = strel('disk',8); %expansion 1 parameters
se2 = strel('disk',3); %expansion 2 parameters
thresh = 1000; %pixel size threshold for objects too small

%define expansion for peripheral zone
periph_distance = 100; % distance in pixels
se_periph=strel('disk', periph_distance);

for p=1:length(points)
    
    % load data
    point=points(p);
    disp(['point',num2str(point)]);
    raw_mask = imread([path,'/Point',num2str(point),'/SMAmask.tiff']);
    figure; imagesc(raw_mask); title('Raw Mask'); plotbrowser on;

    %connect objects and remove objects below threshold
    connected_mask=imclose(raw_mask,se1); %modify se1 for more or less expansion
    figure; imagesc(connected_mask); title('Expanded Mask'); plotbrowser on;
    filtered_mask=MIBIremoveObjectsfromMask(connected_mask,thresh); %modify thresh to tune object sizes
    figure; imagesc(filtered_mask); title('Filtered Mask'); plotbrowser on;
    
    %fill in holes in connected objects
    filled_mask=bwdist(filtered_mask)<=5;
    filled_mask=imfill(filled_mask,'holes');
    figure; imagesc(filled_mask); title('Filled Mask'); plotbrowser on;
    
    % dilate mask to create peripheral zone mask
    expanded_mask=imdilate(filled_mask,se_periph);
    periph_mask=expanded_mask-filled_mask;
    figure; imagesc(periph_mask); title('Peripheral Mask'); plotbrowser on;
    
    % save modified and original masks as matlab variable directory
    imwrite(filled_mask,[path,'/Point',num2str(point),'/region_mask.tif'],'tif');
    imwrite(periph_mask,[path,'/Point',num2str(point),'/peripheral_mask.tif'],'tif');
    save([path,'/Point',num2str(point),'/regionmasks.mat'],'filled_mask','periph_mask');
end
        

