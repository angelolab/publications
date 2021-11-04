function [massDS,countsAllS,countsAllSFilt,baselineVal,baselineValFilt,peaks,peaksFiltered,ratioAll,ratioFiltered,ratioFilteredToNonFiltered,depthsRemovedPerPixel,totalIon,totalIonFilt]= MibiAnalizeSpectra(processedDataDir,pointNumber,depthStart,depthEnd,makePlots,depthProfile,chemicalImage)
% function [countsAllS,countsAllSFilt,baselineLoc,baselineVal,baselineLocFilt,baselineValFilt,peaks,peaksFiltered,ratioAll,ratioFiltered,ratioFilteredToNonFiltered]= MibiAnalizeSpectra(pointNumber,depthStart,depthEnd,makePlots,channelsToPlot,width2baseParam,removeChannels,pathSt)
% function loads spectral data, finds peaks and baseline and plots these                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               as a
% function of plane

if depthProfile
    pathSt = ['Point',num2str(pointNumber),'/RowNumber0/Depth_Profile0/Depth',num2str(depthStart-1)];
elseif chemicalImage
    pathSt = ['Point',num2str(pointNumber),'/RowNumber0/Chemical_Image0'];
else
    error ('If processing msdf, either depth profile or chemical image should be specified');
end
tmp = load([pathSt,'/counts.mat']);
massDS = tmp.massDS;
peaks = zeros(length(massDS),depthEnd);
peaksFiltered = zeros(length(massDS),depthEnd);
baselineVal = zeros(length(massDS),depthEnd);
baselineValFilt = zeros(length(massDS),depthEnd);

countsAllS = cell(1,depthEnd);
countsAllSFilt = cell(1,depthEnd);
baselineAllS = cell(1,depthEnd);
baselineAllSFilt = cell(1,depthEnd);
baselineVal = nan(length(massDS),depthEnd);
baselineValFilt = nan(length(massDS),depthEnd);
depthsRemovedPerPixel= zeros(size(tmp.counts,1));
totalIon=zeros(size(tmp.counts,1),size(tmp.counts,2),depthEnd);
totalIonFilt=zeros(size(tmp.counts,1),size(tmp.counts,2),depthEnd);

for d=depthStart:depthEnd
    if depthProfile
        pathSt = ['Point',num2str(pointNumber),'/RowNumber0/Depth_Profile0/Depth',num2str(d-1)];
    end
    d
    tmp = load([pathSt,'/counts.mat']);
    countsAllS{d} = tmp.counts;
    countsAllSFilt{d} = tmp.countsFilt;
    baselineAllS{d} = tmp.baseline;
    baselineAllSFilt{d} = tmp.baselineFilt;
    peaks(:,d) = tmp.peaks;
    peaksFiltered(:,d) = tmp.peaksFilt;
    baselineVal(:,d) = tmp.baselineVal;
    baselineValFilt(:,d) = tmp.baselineValFilt;
    totalIon(:,:,d) = tmp.totalIon;
    totalIonFilt(:,:,d) = tmp.totalIonFilt;
end

baseWidth = massDS.BaselineEnd - massDS.BaselineStart;
peakWidth = massDS.Stop - massDS.Start;
width2baseParam = peakWidth./baseWidth;

baselineVal = baselineVal.*width2baseParam;
baselineValFilt = baselineValFilt.*width2baseParam;

% get ratio of peaks and baseline
ratioAll = peaks./baselineVal;
ratioFiltered = peaksFiltered./baselineValFilt;
ratioFilteredToNonFiltered = log2(ratioFiltered./ratioAll);

% save
%save (['Point',num2str(pointNumber),'/ParamsSpectra.mat'],'totalIonFilt','totalIon','massDS','countsAllS','countsAllSFilt','baselineLoc','baselineVal','baselineLocFilt','baselineValFilt','peaks','peaksFiltered','ratioAll','ratioFiltered','ratioFilteredToNonFiltered','depthsRemovedPerPixel','-v7.3');
save ([processedDataDir,'/Point',num2str(pointNumber),'/ParamsSpectra.mat'],'totalIonFilt','totalIon','massDS','countsAllSFilt','baselineValFilt','peaks','peaksFiltered','depthsRemovedPerPixel','-v7.3');

%% plots
if makePlots
    disp('Plotting baseline and peaks data');
    
%    % DEBUG: plot total ion intensity with depth
%     f=figure;
%     for d=depthStart:depthEnd
%         subplot(5,ceil(depthEnd/5),d);
%         imagesc(totalIonFilt(:,:,d));
%     end
%     saveas(f,['Point',num2str(pointNumber),'/Results/AllDepths_totalIonFiltered.fig']);
%     close(f);
%         
%     % plot selected channels of non-filtered data
%     [tf,channelsLoc] = ismember(channelsToPlot,massDS.Label);
%     for i=1:length(channelsToPlot)
%         f=figure;
%         for d=depthStart:depthEnd
%             subplot(5,ceil(depthEnd/5),d);
%             imagesc(countsAllS{d}(:,:,channelsLoc(i)));
%         end
%         saveas(f,['Point',num2str(pointNumber),'/Results/AllDepths_',channelsToPlot{i},'.fig']);
%         close(f);
%     end

    % plot baseline value and location for ll markers for non-filtered data
    f0=figure;
    plot(baselineVal');
    xlabel('planes');
    ylabel('baseline value');
    legend(massDS.Label);
    grid on;
    saveas(f0,[processedDataDir,'/Point',num2str(pointNumber),'/Results/shiftsInBaseLineAlongPlanes.fig']);

%     % DEBUG: plot shift in baseline position       
%     f1=figure;
%     plot(baselineLoc');
%     xlabel('planes');
%     ylabel('baseline location');
%     legend(massDS.Label);
%     grid on;
%     saveas(f1,['Point',num2str(pointNumber),'/Results/shiftsInBaseLinePosAlongPlanes.fig']);

    % plot baseline value and location for all markers for filtered data
    f2=figure;
    plot(baselineValFilt');
    xlabel('planes');
    ylabel('baseline value');
    title('Filtered Image');
    legend(massDS.Label);
    grid on;
    saveas(f2,[processedDataDir,'/Point',num2str(pointNumber),'/Results/shiftsInBaseLineAlongPlanesForFilteredImage.fig']);

%     % DEBUG: plot shift in baseline position      
%     f3=figure;
%     plot(baselineLoc');
%     xlabel('planes');
%     ylabel('baseline location');
%     title('Filtered Image');
%     legend(massDS.Label);
%     grid on;
%     saveas(f3,['Point',num2str(pointNumber),'/Results/shiftsInBaseLinePosAlongPlanesForFilteredImage.fig']);

    % plot counts
    f4=figure;
    plot(peaks');
    xlabel('planes');
    ylabel('peak data');
    title('Non-Filtered Image');
    legend(massDS.Label);
    grid on;
    saveas(f4,[processedDataDir,'/Point',num2str(pointNumber),'/Results/peaksAlongPlanes.fig']);

    % plot counts
    f5=figure;
    plot(peaksFiltered');
    xlabel('planes');
    ylabel('peak data');
    title('Filtered Image');
    legend(massDS.Label);
    grid on;
    saveas(f5,[processedDataDir,'/Point',num2str(pointNumber),'/Results/peaksAlongPlanesFilteredData.fig']);

    % DEBUG: plot ratio of peak to baseline
    f6=figure;
    plot(ratioAll');
    xlabel('planes');
    ylabel('ratio of peak to baseline');
    title('Non-Filtered Image');
    legend(massDS.Label);
    grid on;
    saveas(f6,[processedDataDir,'/Point',num2str(pointNumber),'/Results/ratioOfPeaksToBaselineAlongPlanes.fig']);

    f7=figure;
    plot(ratioFiltered');
    xlabel('planes');
    ylabel('ratio of peak to baseline');
    title('Filtered Image');
    legend(massDS.Label);
    grid on;
    saveas(f7,[processedDataDir,'/Point',num2str(pointNumber),'/Results/ratioOfPeaksToBaselineAlongPlanesForFilteredData.fig']);

    % plot one over the ratio, to get the paercentage of positive signal that
    % needs to be removed
    f8=figure;
    plot(1./ratioFiltered');
    xlabel('planes');
    ylabel('ratio of baseline to peak');
    title('Filtered Image');
    legend(massDS.Label);
    grid on;
    saveas(f8,[processedDataDir,'/Point',num2str(pointNumber),'/Results/ratioOfBaselineToPeaksAlongPlanesForFilteredData.fig']);


    % compare ratios of filtered and non-filtered data. Do we improve our
    % signal to noise?
    f9=figure;
    plot(ratioFilteredToNonFiltered');
    xlabel('planes');
    ylabel('Log2(Signal-to-noise in filtered data / Signal-to-noise in non-filtered data)');
    title('Comparing filtered and non-filtered data');
    legend(massDS.Label);
    grid on;
    line([1,depthEnd],[0,0],'LineWidth',2,'Color',[0 0 0]);
    saveas(f9,[processedDataDir,'/Point',num2str(pointNumber),'/Results/ratioOfSignalToNoiseInFilteredAndNonfilteredData.fig']);

    % plot estimate of bad counts in non-filtered data and the amount to
    % remove according to the bg channels
    goodEstAll = peaks.*(1- (1./ratioAll));
    badEstAll = peaks.*(1./ratioAll);
    removedInFilter = peaks-peaksFiltered;
    removedInFilterPercent = (peaks-peaksFiltered)./peaks;
    f10=figure;
    plot(badEstAll');
    hold on;
    plot(removedInFilter')
    xlabel('planes');
    ylabel('estimate bad counts non filtered and how much to remove');
    legend([massDS.Label;massDS.Label]);
    grid on;
    saveas(f10,[processedDataDir,'/Point',num2str(pointNumber),'/Results/estimatedCountsToFilterAndAmountFilteredByGold.fig']);

    
    % same, but with percentages
    f11=figure;
    plot(1./ratioAll');
    hold on;
    plot(removedInFilterPercent')
    xlabel('planes');
    ylabel('estimate baseline non filtered and percentage to remove by gold');
    legend([massDS.Label;massDS.Label]);
    grid on;
    saveas(f11,[processedDataDir,'/Point',num2str(pointNumber),'/Results/baselineAndPercentFilteredByGold.fig']);

    
end