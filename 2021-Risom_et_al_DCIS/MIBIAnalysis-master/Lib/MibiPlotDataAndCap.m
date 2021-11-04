function MibiPlotDataAndCap(data,cap,titlestr)
% function MibiPlotDataAndCap(data,cap)
% function plots the data and sets any value larger than cap to cap

currdata = data;
currdata(currdata>cap) = cap;
figure;
imagesc(currdata);
title(titlestr);