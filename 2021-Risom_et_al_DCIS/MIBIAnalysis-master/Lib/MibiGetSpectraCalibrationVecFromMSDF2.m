function specificSpectraVec = MibiGetSpectraCalibrationVecFromMSDF2(fileMSDFAccum,paramS,massDS,spectraVec)
% function spectraVec = MibiGetSpectraCalibrationVecFromMSDF2(fileMSDFAccum,paramS,massDS,spectraVec)
% function gets a spectra vec listing a correspondance between the general
% region of two defined tof peaks and their masses. It reads an msdf
% spectrum and returns the exact values of the spectra peaks in this msdf
% file.

tofData = MibiParseMsdfSpectrum(fileMSDFAccum,paramS,massDS);

% bin data
edges = [0:paramS.TimeResolution:10000];
specHist = histcounts(tofData,edges);

% find the highest peak within 30 bins to the value that we inserted
% (15 above and 15 below)
n=15;
% find the bin of the inserted value
[~,givenValLoc1] = min(abs(edges-spectraVec(1)));
[~,givenValLoc2] = min(abs(edges-spectraVec(3)));

neighborhood1 = specHist([givenValLoc1-n:givenValLoc1+n]);
edgesNeighborhood1 = edges([givenValLoc1-n:givenValLoc1+n]);
neighborhood2 = specHist([givenValLoc2-n:givenValLoc2+n]);
edgesNeighborhood2 = edges([givenValLoc2-n:givenValLoc2+n]);

[maxVal1,indMaxLoc1] = max(neighborhood1);
[maxVal2,indMaxLoc2] = max(neighborhood2);

specificSpectraVec = spectraVec;
specificSpectraVec(1) = edgesNeighborhood1(indMaxLoc1);
specificSpectraVec(3) = edgesNeighborhood2(indMaxLoc2);


