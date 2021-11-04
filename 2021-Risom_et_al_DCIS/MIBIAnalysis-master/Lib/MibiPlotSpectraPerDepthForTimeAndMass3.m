function MibiPlotSpectraPerDepthForTimeAndMass3(fileNameXML,fileNameMass,pointNumber,depthStart,depthEnd,figureCompareTOFFileName,figureCompareMassFileName,calibrateSpectra,calibrateSpectraPerDepth,spectraVec)
% function MibiPlotSpectraPerDepth
% plots overlapping spectra for specified depths to see if shifts occured during the run

% get mass calibration parameters
paramS = MibiReadRunXml(fileNameXML,pointNumber);

% get mass edges
massDS = MibiReadMassData(fileNameMass);
massNum = size(massDS,1);

% generate labels
histLabels = cell(depthEnd-depthStart+1,1);
for d=depthStart:depthEnd
    histLabels{d} = ['depth',num2str(d)];
end

% read data
tofPerDepth = cell(depthEnd-depthStart+1,1);
massPerDepth = cell(depthEnd-depthStart+1,1);

for d=depthStart:depthEnd
    disp(['processing depth',num2str(d-1)]);
    pathSt = ['Point',num2str(pointNumber),'/RowNumber0/Depth_Profile0/Depth',num2str(d-1)];
    fileMSDF=[pathSt,'/ImageAccum.msdf'];
    tofPerDepth{d} = MibiParseMsdfSpectrum(fileMSDF,paramS,massDS);
    
    % bin data
    edges = [0:paramS.TimeResolution:20000];
    specHist = histcounts(tofPerDepth{d},edges);
    
    % find the highest peak within 30 bins to the value that we inserted
    % (15 above and 15 below)
    n=15;
    % find the bin of the inserted value
    [~,givenValLoc1] = min(abs(edges-spectraVec(1)));
    [~,givenValLoc2] = min(abs(edges-spectraVec(3)));
    
    neighborhood1 = specHist([givenValLoc1-n:givenValLoc1+n]);
    edgesNeighborhood1 = edges([givenValLoc1-n:givenValLoc1+n]);
    neighborhood2 = specHist([givenValLoc2-n:givenValLoc2+n]);
    edgesNeighborhood2 = edges([givenValLoc2-n:givenValLoc2+n]);
    
    [maxVal1,indMaxLoc1] = max(neighborhood1);
    [maxVal2,indMaxLoc2] = max(neighborhood2);
 
%     %% get spectra my max
%     specificSpectraVec = spectraVec;
%     specificSpectraVec(1) = edgesNeighborhood1(indMaxLoc1);
%     specificSpectraVec(3) = edgesNeighborhood2(indMaxLoc2);

    
    %% get spectra by mean
    % get the index for the mean of the peak (considered as 50% of max and
    % above)
    maxPos1 = edgesNeighborhood1(indMaxLoc1);
    maxPos2 = edgesNeighborhood2(indMaxLoc2);
    
    % find the bin of the max value
    [~,indMax1] = min(abs(edges-maxPos1));
    [~,indMax2] = min(abs(edges-maxPos2));
    
    max1NeighborsInds = ([indMax1-n:indMax1+n]);
    max2NeighborsInds = ([indMax2-n:indMax2+n]);
    max1NeighborsVals = specHist(max1NeighborsInds);
    max2NeighborsVals = specHist(max2NeighborsInds);
    
    % find the vals above 50% of max
    valsAbove50PercentAroundMax1 = max1NeighborsVals>0.5*maxVal1;
    meanPosAbove50PercentAroundMax1 = mean(edges(max1NeighborsInds(valsAbove50PercentAroundMax1)));
    valsAbove50PercentAroundMax2 = max2NeighborsVals>0.5*maxVal2;
    meanPosAbove50PercentAroundMax2 = mean(edges(max2NeighborsInds(valsAbove50PercentAroundMax2)));
    
    specificSpectraVec = spectraVec;
    specificSpectraVec(1) = meanPosAbove50PercentAroundMax1;
    specificSpectraVec(3) = meanPosAbove50PercentAroundMax2;
    
    %% new spectra
    if (calibrateSpectra == 1)
        sol= MibiCalibrateSpectrum (specificSpectraVec);
        paramS.MassGain=eval(sol.a);
        paramS.MassOffset=eval(sol.b);
        if length(paramS.MassGain)>1
            paramS.MassGain=paramS.MassGain(1);
            paramS.MassOffset=paramS.MassOffset(1);
        end
    end
    massPerDepth{d} = (( tofPerDepth{d} / 1000 - paramS.MassOffset ) / paramS.MassGain).^2;
end

% bin data tof
edges = [0:paramS.TimeResolution:20000];
dataVec = zeros(length(edges)-1,depthEnd-depthStart+1);
for d=depthStart:depthEnd
    dataVec(:,d) = histcounts(tofPerDepth{d},edges);
end

% plot as function of time
h=figure;
c=colormap(jet(depthEnd-depthStart+1));
for d=depthStart:depthEnd
    hold on;
    plot(edges(1:end-1),dataVec(:,d),'Color', c(d,:));
end
legend(histLabels);
grid on;

saveas(h,figureCompareTOFFileName);

% bin data mass
edgesMass = [0:0.01:200];
dataVecMass = zeros(length(edgesMass)-1,depthEnd-depthStart+1);
for d=depthStart:depthEnd
    dataVecMass(:,d) = histcounts(massPerDepth{d},edgesMass);
end

%% checking the mass calibration
% for each expected peak get the max
measuredPeaks = zeros(length(massDS),depthEnd-depthStart+1);
maxPeaks = zeros(length(massDS),depthEnd-depthStart+1);
for d=depthStart:depthEnd
    for m=1:length(massDS)
        [tfStart, locStart]=ismember(massDS.Start(m),round(edgesMass,2));
        [tfEnd, locEnd]=ismember(massDS.Stop(m),round(edgesMass,2));
        currPeakEdges = edgesMass([locStart:locEnd]);
        currPeakVals = dataVecMass([locStart:locEnd],d);
        [currMax,currMaxInd] = max(currPeakVals);
        maxPeaks(m,d) = currPeakEdges(currMaxInd);
        above50PercentLog = currPeakVals>0.5*currMax;
        currPeakEdgesAbove50 = currPeakEdges(above50PercentLog);
        measuredPeaks(m,d) = mean(currPeakEdgesAbove50);
    end
end

% plot as a function of mass
f=figure;
c=colormap(jet(depthEnd-depthStart+1));
for d=depthStart:depthEnd
    hold on;
    plot(edgesMass(1:end-1),dataVecMass(:,d),'Color', c(d,:));
end
legend(histLabels);
grid on;

% plot without zero vals
f=figure;
c=colormap(jet(depthEnd-depthStart+1));
for d=depthStart:depthEnd
    hold on;
    dataVecTrunc=dataVecMass(:,d);
    plotInds = dataVecTrunc>0;
    plot(edgesMass(plotInds),dataVecTrunc(plotInds),'Color', c(d,:),'LineWidth',2);
end
xlabel('mass');
ylabel('Counts');

saveas(f,figureCompareMassFileName);
