function  [sol] = MibiCalibrateSpectrum (spectraVec)
% function MibiCalibrateSpectrum (spectraVec)
% function receives two time values and their corresponding masses and
% computes the time resolution (a) and offset (b) to calibrate the spectrum
% according to: m(daltons) = a*t^2 + b

syms a b c

t1= spectraVec(1);
m1= spectraVec(2);
t2= spectraVec(3);
m2= spectraVec(4);

sol = solve(((t1/1000 - b)/a)^2 - m1 , ((t2/1000 - b)/a)^2 - m2);