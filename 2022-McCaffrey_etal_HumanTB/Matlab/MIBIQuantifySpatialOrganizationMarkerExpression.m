% % MIBIQuantifySpatialOrganizationMarkerExpression.m
% Author: Leeat Keren, adapted to TB cohort by Erin McCaffrey

% Get enrichment of cells positive for certain proteins to sit together or
% not. Quantify using a defined pixel distance. For each sample get all the
% interactions (defined as pixel distance) between cells positive for protein X
% and protein Y. The use the bootstrapping approach to permute the labels randomly

 
%% 1. Define paths and read in data

path = [cd,'/data/']; % path to sc data and marker thresholds
pathSegment = [cd,'/data/segmentation/']; % path to distance matrices
% TB points to analyze
points = [21,84,42,88,28,89,90,91,94,95,96,97,14,15,98,99,6,7,33,34,...
          26,27,40,61,47,48,54,55,92,93]; 
      
% import data
sc_data = [path,'allTB-sarcoid-scdata_matlab.csv'];
opts = detectImportOptions(sc_data);
dataAll = readtable(sc_data,opts);

%converting to matrix for easy indexing
dataAllMat=table2cell(dataAll);
dataAllMat=cell2mat(dataAllMat(:,1:44)); %exclude str annotation

%% 2. Define marker pairs

markerInds = [6:42]; %indices of the markers in matrix
dataMarkers = dataAllMat(:,markerInds); %subset of data matrix including only markers
markerTitles = dataAll(:,markerInds).Properties.VariableNames; %names of the markers
markerNum = length(markerTitles); %number of markers to compare

%% 3. Set up permutation conditions

BootstrapNum = 1000; %number of permutations for bootsrapping to create null
distLim = 100; %cutoff to define the cells that are interacting/close in pixels
closeNum = zeros(markerNum,markerNum); %matrix for number of interactions for each marker pair
closeNumRand = zeros(markerNum,markerNum,BootstrapNum); %matrix for null number of interactions

%% 4. Define some useful variables for indexing

patientIdx = 1; %column with patient/point label
cellLabelIdx = 2; %column with cell label

%% 5. Read in marker thresholds (custom values in same order as matrix columns)

markerThresh=readtable([path,'markerThresholds.csv']); 
thresh_vec = markerThresh.Threshold;

%% 6. Run enrichment
for i=1:length(points)
    point=points(i);
    % load relevant data
    disp(['Working on point:',num2str(point)]);
    load([pathSegment,'/Point',num2str(point),'_cellDistances.mat']);
    % get data for current patient
    patientInds = dataAllMat(:,patientIdx) == point; %get row indices of current point
    patientData = dataAllMat(patientInds,:); %get just the data for the current point
    patientDataMarkers = dataMarkers(patientInds,:); %get just the marker data for the current point
    % go over markers
    for j=1:markerNum
        marker1_thresh = thresh_vec(j); %get positiviy threshold for marker
        marker1PosInds = (patientDataMarkers(:,j) > marker1_thresh ); %indices of cells positive for that marker
        marker1PosLabels = patientData(marker1PosInds,cellLabelIdx); %labels of cells positive for that marker
        marker1Num = length(marker1PosLabels); %number of cells positive for that marker
        % iterate over all other markers + curr marker
        parfor k=1:markerNum
            marker2_thresh = thresh_vec(k);
            marker2PosInds = ( patientDataMarkers(:,k) > marker2_thresh ); %indices of cells positive for that marker
            marker2PosLabels = patientData(marker2PosInds,cellLabelIdx); %labels of cells positive for that marker
            marker2Num = length(marker2PosLabels); %number of cells positive for that marker
            truncDistMat = distancesMat(marker1PosLabels,marker2PosLabels); %create a distance matrix of just pos v pos cells
            % turn to binary
            truncDistMatBin = zeros(size(truncDistMat));
            truncDistMatBin(truncDistMat<distLim) = 1; %convert to binary of interacting or not
            % record interaction num
            closeNum(j,k) = sum(sum(truncDistMatBin))  %add number of interactions to matrix of 'real' interactions
            % randomize to get null distribution
            for r=1:BootstrapNum
                % get random labels for marker 1
                marker1LabelsRand = datasample(patientData(:,cellLabelIdx),marker1Num);
                % get random labels for marker 2
                marker2LabelsRand = datasample(patientData(:,cellLabelIdx),marker2Num);
                randTruncDistMat = distancesMat(marker1LabelsRand,marker2LabelsRand); % distance matrix of the random labels
                % turn to binary
                randTruncDistMatBin = zeros(size(randTruncDistMat)); %convert to binary of interacting or not
                randTruncDistMatBin(randTruncDistMat<distLim) = 1; 
                % record interaction num
                closeNumRand(j,k,r) = sum(sum(randTruncDistMatBin)); %add number of interactions to matrix of 'null' interactions
            end
        end
    end

    % Create a vector for the z-scoring (and needed params mean + std)
    z = zeros(markerNum); %z-score
    muhat = zeros(markerNum); %mean
    sigmahat = zeros(markerNum); %std deviation
    
    % Create a vector for the empirical p-value 
    p = zeros(markerNum,markerNum,2);

    % Get the enrichment z-score for each marker
    for j=1:markerNum
        for k=1:markerNum
            tmp= reshape(closeNumRand(j,k,:),BootstrapNum,1);
            [muhat(j,k),sigmahat(j,k)] = normfit(tmp);
            z(j,k) = (closeNum(j,k)-muhat(j,k))/sigmahat(j,k);
            p(j,k,1) = (1+(sum(tmp>=closeNum(j,k))))/(BootstrapNum+1);
            p(j,k,2) = (1+(sum(tmp<=closeNum(j,k))))/(BootstrapNum+1);
        end
    end
    
    % adjust p-values using FDR 0.05 (Inf or NaN will have p value 1)
    p_summary = p(:,:,1);
    for j=1:markerNum
        for k=1:markerNum
            % if interaction is enriched +z grab first p-value
            if (z(j,k) > 0)
                p_summary(j,k) = p(j,k,1);
            % if interaction is avoided -z grab second p-value
            elseif (z(j,k) < 0) 
                p_summary(j,k) = p(j,k,2);
            end
        end
    end
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p_summary,.05,'pdep','yes');

    % save results
    save([pathSegment,'/Point',num2str(point),'_spatialAnalysis.mat'],...
        'closeNum','closeNumRand','z','muhat','sigmahat','p','adj_p','h','markerTitles');

end


%% 7. Optionally visualize individual heatmaps
for i=1:length(points)
    point=points(i);
    disp(['Working on point:',num2str(point)]);
    load([pathSegment,'/Point',num2str(point),'_spatialAnalysis.mat'])
    % plot
    zplot = z;
    zplot(isnan(zplot)) = 0;
    zplot(isinf(zplot)) = 0;
    hmap=clustergram(zplot,'RowLabels', markerTitles,'ColumnLabels',...
        markerTitles,'Colormap', 'redbluecmap','DisplayRange', 20,'DisplayRatio', 0.2);
    addTitle(hmap, ['Point ',num2str(point)])
    fig=plot(hmap);
end

