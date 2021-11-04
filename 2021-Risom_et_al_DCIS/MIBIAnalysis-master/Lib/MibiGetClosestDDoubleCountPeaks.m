function [IDX, closestDvecBA] = MibiGetClosestDDoubleCountPeaks(markerApeaks,markerBpeaks,K)
% function MibiGetClosestDDoubleCountPeaks(dataA,dataB,K)
% function calculates for dataA and dataB the distance to the K closest
% neighbours. Peaks with value above 1 get counted several times

[xA yA]=find(markerApeaks);
[xB yB]=find(markerBpeaks);

% duplicate intances of positive pixels in markerA according to their peak
% intensity
xAexpand = repelem(xA,markerApeaks(markerApeaks>0));
yAexpand = repelem(yA,markerApeaks(markerApeaks>0));

[IDX,closestDvecBA] = knnsearch([xAexpand, yAexpand],[xB yB],'K',K);