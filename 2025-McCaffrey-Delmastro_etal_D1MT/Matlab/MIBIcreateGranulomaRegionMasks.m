% MIBIcreateGranulomaRegionMasks.m
% Date created: 190620
% Author: Erin McCaffrey 

% MIBI automatically assign pixels as being in a myeloid-rich region or the
% peripheral region.

%define paths and points to analyse
path = '/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/Reviewer_Experiments/D1MT-cohort/Cohort/single-roi-master/no_noise/';
massDS = MibiReadMassData('/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/Reviewer_Experiments/D1MT-cohort/Cohort/single-roi-master/GranulomaPanel.csv');
points = [8,9,13,14,15,16,18,25,27,29,30,31,39,40,41,1,2,10,11,12,19,...
    20,21,22,23,24,28,32,33,34,35,36,37,38]; %cohort data
imSize=1024;

%select channels from spatiel enrichment to use for analysis
myeloid_channels={'IDO'};

%get indices of channels
[~,myeloidInds]=ismember(myeloid_channels,massDS.Label);

%define masking parameters
t=0.1;
gausRad=5;
cap=10;
se1 = strel('disk',10); %expansion for myeloid
se2 = strel('disk',3); %expansion for myeloid
thresh = 10000;
periph_distance = 100; % pixel distance to expand for 'peripheral zone

for p=1:length(points)
    % load data
    point=points(p);
    disp(['point',num2str(point)]);
    load([path,'/Point',num2str(point),'/dataNoAgg.mat'])
    
    %create myeloid mask
    combined=MIBIcreateSummedChannel(myeloidInds,countsNoNoiseNoAgg);
    figure; imagesc(combined); title('Myeloid Data'); plotbrowser on;
    myeloid_mask=MibiGetMask(combined,cap,t,gausRad);
    figure; imagesc(myeloid_mask); title('Unprocessed Myeloid Mask'); plotbrowser on;
    %connect objects and remove objects below threshold
    myeloid_mask=imclose(myeloid_mask,se1);
    figure; imagesc(myeloid_mask); title('Expanded Myeloid Mask'); plotbrowser on;
    myeloid_mask=MIBIremoveObjectsfromMask(myeloid_mask,thresh);
    figure; imagesc(myeloid_mask); title('Filtered Myeloid Mask'); plotbrowser on;
    %fill in holes in connected objects
    myeloid_mask=bwdist(myeloid_mask)<=2;
    myeloid_mask=imfill(myeloid_mask,'holes');
    figure; imagesc(myeloid_mask); title('Complete Myeloid Mask'); plotbrowser on;
    % dilate mask before active contouring
    expanded_mask=imdilate(myeloid_mask,se2);
    figure; imagesc(expanded_mask); title('Expanded Myeloid Mask'); plotbrowser on;
    final_mask=imfill(activecontour(combined,expanded_mask,200, 'edge','SmoothFactor',10),'holes');
    figure; imagesc(final_mask); title('Final Mask'); plotbrowser on;
    close all;
    imwrite(uint16(final_mask),[path,'/Point',num2str(point),'/myeloid_mask.tif'],'tif');
end
        

