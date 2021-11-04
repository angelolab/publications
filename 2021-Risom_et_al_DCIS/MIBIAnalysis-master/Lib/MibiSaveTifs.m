function MibiSaveTifs (savePath, dataMat, channelNames)
% function MibiSaveTifs (folder, dataMat, chanelNames)
% function gets 3d Mibi data and saves it as individual tif files

if 7~=exist(savePath,'dir')
    mkdir(savePath);
end

for i=1:length(channelNames)
    data = uint16(dataMat(:,:,i));
    filename = [savePath,filesep,channelNames{i},'.tiff'];
    imwrite(data,filename);
    t = Tiff(filename, 'r+');
    setTag(t, 'PageName', channelNames{i});
    rewriteDirectory(t);
    close(t);
end