d = laserdata();
Jmin = min([d.exp.JS1x; d.exp.JS2x]);
Jmax = max([d.exp.JS1x; d.exp.JS2x]);
JArr = linspace(Jmin, Jmax, 50).';

% d.sig = 20e-3;
% x0  = [0.1190    0.0414    5.3972    1.1767];
% d.eta = x0(1);
% d.taun = x0(2);
% d.Z = x0(3);
% d.E0J0 = x0(4);

d.sig = 20e-3;
x0  = [0.1058    0.0897    7.8422    1.1788];
d.eta = x0(1);
d.taun = x0(2);
d.Z = x0(3);
d.E0J0 = x0(4);

%save spectra to files at each injection current
cfg.writeSpectra = 1;

%integration parameters
%range of QD energies in standard deviation units (sigma)
%   1 accounts for 68.7% of all QDs
%   2 - 95.45%
%   3 - 99.73%
%however, since most of the action happens on the very edge
%   of QD energy distribution
cfg.rangeE = 6;
%number of distretization points
cfg.ptsE = 300;

%how long laser is expected to reach equilibrium, nsec
cfg.integT = 1e3; 


out = calc(JArr, d, cfg);
[s, powerFine] = powerScale(out, d);
fprintf('Power scaling fine %f\n', powerFine)
out.S1Arr = out.S1Arr * s;
out.S2Arr = out.S2Arr * s;

csvwrite('data/state.csv', [out.JArr/1e3, out.efArr, out.nArr, out.S1Arr, ...
    out.S2Arr, out.Gain1, out.Gain2, out.e0Arr, out.e1Arr, out.e2Arr, ...
    out.GmaxE, out.GainFWHM])

draw(out, d)