function countsNoBg = MibiRemoveBackgroundByMaskAllChannels(countsAllSFiltCRSum,mask, removeVal)
% function MibiRemoveBackgroundByMaskAllChannels
% function receives a matrix of counts and a logical 2d-mask, and reduces
% values where mask is positive. The amount of reduction is set by
% removeVal. If it is not provided, the function will zero these regions

channelNum = size(countsAllSFiltCRSum,3);
mask3d = repmat(mask,1,1,channelNum);
countsNoBg = countsAllSFiltCRSum;
if ~exist('removeVal')
    countsNoBg(mask3d) = 0;
else
    countsNoBg(mask3d) = countsNoBg(mask3d)-removeVal;
    countsNoBg(countsNoBg<0) = 0;
end