function [intND] = MibiGetIntNormDist(dataA,dataB,K,n2joinStart,n2joinEnd)
% function MibiGetIntNormDist(dataA,dataB,K)
% function calculates for dataA and dataB the distance to the K closest
% neighbours. Peaks with value above 1 get counted several times according
% to their peak value. A score vector of the Intensity-normalized-distance
% is returned. The score is calculated by the average distance to the
% neighbours n2joinStart to n2joinEnd
% For example: MibiGetIntNormDist(dataA,dataB,5,2,5) will get the distances
% to the 5 nearest neighbours, counting peaks>1 more than once. It will
% then return the average distance to neighbours 2-5.

[IDX, closestDvecBA] = MibiGetClosestDDoubleCountPeaks(dataA,dataB,K);
if isempty(closestDvecBA)
    intND=[];
else
    intND = mean(closestDvecBA(:,[n2joinStart:n2joinEnd]),2);
end
