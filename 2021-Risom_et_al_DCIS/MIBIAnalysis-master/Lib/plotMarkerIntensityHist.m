function f=plotMarkerIntensityHist (allMarkers, ImageNames)
% function plotMarkerIntensityHist (allMarkers, ImageNames)
% function plots histograms of intensity for all markers

markerNum=length(ImageNames);
colNum=5;

f=figure;
for i=1:markerNum
    subplot(ceil(markerNum/colNum),colNum,i);
    currData=allMarkers(:,:,i);
    currData= currData(:);
    currData= currData(currData>0);
    [binY,binX]=hist(currData(:));
    bar(binX,log10(binY));
    %bar(binX,binY);
    title(ImageNames(i));
    axis tight;
end
    
