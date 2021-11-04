function filteredData = MibiFilterImageByNNThreshold(data,NNVec,nnValToRemov)
% function filteredData = MibiFilterImageByNNThreshold(data,NNVec,nnValToRemov)
% function receives an image and filters it according to the nearest
% neighbour results

filteredData = data;
posSignal=find(filteredData);
if ~isempty(posSignal)
    [nnScorS, inds] = sort(NNVec,'descend');
    pos=nnScorS>nnValToRemov;
    indsToRemove = inds(pos);
    pos2Remove=posSignal(indsToRemove);
    filteredData(pos2Remove) = 0;
end
