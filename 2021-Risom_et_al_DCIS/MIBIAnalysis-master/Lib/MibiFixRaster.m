function [countsAllSFiltNew,totalIonNew,totalIonFiltNew]= MibiFixRaster(countsAllSFilt,totalIon,totalIonFilt, fixRasterCols, fixRasterRows)
% function fixes cases in which raster scanner goes crazy

% make a matrix out of the carbon channel
fixChannel = 5;
channelMat = zeros(size(totalIon));
for i= 1:size(totalIon,3)
    channelMat(:,:,i) = countsAllSFilt{i}(:,:,fixChannel);
end
channelMat(isnan(channelMat)) = 0;

if fixRasterCols
    totalIonDeriv= channelMat(:,[2:end],:)-channelMat(:,[1:end-1],:);
    sumDeriv = nansum(abs(totalIonDeriv),1);
    sumDeriv2 = reshape(sumDeriv,size(totalIonDeriv,2),size(totalIonDeriv,3));

    secondDeriv = sumDeriv2([2:end],:)-sumDeriv2([1:end-1],:);


    % [envHigh, envLow] = envelope(sumDeriv2(:,1),10,'peak');
    % envMean = (envHigh+envLow)/2;
    % 
    % plot([1:size(sumDeriv2,1)],sumDeriv2(:,1), ...
    %      [1:size(sumDeriv2,1)],envHigh, ...
    %      [1:size(sumDeriv2,1)],envMean, ...
    %      [1:size(sumDeriv2,1)],envLow)


    % l=cell(1,size(totalIonDeriv,3));
    % for i=1:size(totalIonDeriv,3)
    %     l{i} = num2str(i);
    % end
    % figure;
    % plot(sumDeriv2);
    % legend(l);

    % get the max values of the derivitive
    [maxDerivVal, maxDerivLoc] = max(secondDeriv);
    % move cols
    countsAllSFiltNew = cell(size(countsAllSFilt));
    totalIonNew = zeros(size(totalIon));
    totalIonFiltNew = zeros(size(totalIonFilt));
    for i=1:length(maxDerivLoc)
        currMat = countsAllSFilt{i};
        newMat = currMat(:,[maxDerivLoc(i)+1:end,1:maxDerivLoc(i)],:);
        countsAllSFiltNew{i} = newMat;
        totalIonNew(:,:,i) = totalIon(:,[maxDerivLoc(i)+1:end,1:maxDerivLoc(i)],i);
        totalIonFiltNew(:,:,i) = totalIonFilt(:,[maxDerivLoc(i)+1:end,1:maxDerivLoc(i)],i);
    end
end

if fixRasterRows
    if fixRasterCols
        % make sure new values are fed into row correction
        fixChannel = 5;
        channelMat = zeros(size(totalIonNew));
        for i= 1:size(totalIon,3)
            channelMat(:,:,i) = countsAllSFiltNew{i}(:,:,fixChannel);
        end
        channelMat(isnan(channelMat)) = 0;
        
        countsAllSFilt = countsAllSFiltNew;
        totalIon = totalIonNew;
        totalIonFilt = totalIonFiltNew;
    end
        
    totalIonDeriv= channelMat([2:end],:, :)-channelMat([1:end-1],:,:);
    sumDeriv = nansum(abs(totalIonDeriv),2);
    %sumDeriv2 = reshape(sumDeriv,size(totalIonDeriv,1),size(totalIonDeriv,3));

    secondDeriv = sumDeriv([2:end],:)-sumDeriv([1:end-1],:);


    % [envHigh, envLow] = envelope(sumDeriv2(:,1),10,'peak');
    % envMean = (envHigh+envLow)/2;
    % 
    % plot([1:size(sumDeriv2,1)],sumDeriv2(:,1), ...
    %      [1:size(sumDeriv2,1)],envHigh, ...
    %      [1:size(sumDeriv2,1)],envMean, ...
    %      [1:size(sumDeriv2,1)],envLow)


    % l=cell(1,size(totalIonDeriv,3));
    % for i=1:size(totalIonDeriv,3)
    %     l{i} = num2str(i);
    % end
    % figure;
    % plot(sumDeriv2);
    % legend(l);

    % get the max values of the derivitive
    [maxDerivVal, maxDerivLoc] = max(secondDeriv);
    % move cols
    countsAllSFiltNew = cell(size(countsAllSFilt));
    totalIonNew = zeros(size(totalIon));
    totalIonFiltNew = zeros(size(totalIonFilt));
    for i=1:length(maxDerivLoc)
        currMat = countsAllSFilt{i};
        newMat = currMat([maxDerivLoc(i)+1:end,1:maxDerivLoc(i)],:,:);
        countsAllSFiltNew{i} = newMat;
        totalIonNew(:,:,i) = totalIon([maxDerivLoc(i)+1:end,1:maxDerivLoc(i)],:,i);
        totalIonFiltNew(:,:,i) = totalIonFilt([maxDerivLoc(i)+1:end,1:maxDerivLoc(i)],:,i);
    end
end

 
