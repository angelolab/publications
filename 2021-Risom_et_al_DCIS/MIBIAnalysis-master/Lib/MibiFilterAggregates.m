function dataNoAgg = MibiFilterAggregates(data, gausRad, t, gausFlag)
% dataNoAgg = MibiRemoveAggregates(data, gausRad, t)
% function receives a 2d counts matrix and removes
% aggregates. Matrix undergoes gaussian kernel with radius gausRad. Aggregates are considered as connected components of size>t
% gausFlag - 0 if no gaussian should be done
% plotFlag - 0 if plots should not be generated

    if ~exist('gausFlag')
        gausFlag = 0;
    end

    dataGaus = data;
    if gausFlag
        dataGaus = imgaussfilt(data,gausRad);
    end
    posdata = imbinarize(dataGaus,0);

    stats = regionprops(posdata,'Area','PixelIdxList');
    dataNoAgg = data;
    for i=1:length(stats)
        if stats(i).Area <t
            dataNoAgg(stats(i).PixelIdxList) = 0;
        end
    end

end