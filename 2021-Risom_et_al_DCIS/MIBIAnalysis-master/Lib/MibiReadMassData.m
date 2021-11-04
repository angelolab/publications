function massDS = MibiReadMassData(fileName)
% function massDS = MibiReadMassData(fileName)
% function reads mass edges and labels from csv-delimited files

massDS = dataset('File',fileName,'Delimiter',',');

