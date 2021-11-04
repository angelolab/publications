function countsNoNoise = MibiFilterAllByNN(countsAllSFiltCRSum,IntNormD,threshVec)
% function filteredData = MibiFilterAllByNN(countsAllSFiltCRSum,IntNormD,threshFile)
% function filters images (removes noise) by NN density. countsAllFiltCRSum
% - cell array of images. IntNormD - cell array of
% intensity-normalized-distances. threshVec - a vector that indicates the
% filtering threshold for each marker in the panel. Values of zero indicate
% that no filtering should be done.

countsNoNoise=zeros(size(countsAllSFiltCRSum));
for i=1:size(countsAllSFiltCRSum,3)
    if threshVec(i)>0
        countsNoNoise(:,:,i) = MibiFilterImageByNNThreshold(countsAllSFiltCRSum(:,:,i),IntNormD{i},threshVec(i));
    else
        countsNoNoise(:,:,i) = countsAllSFiltCRSum(:,:,i);
    end
end