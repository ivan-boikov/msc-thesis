d = laserdata();

%initial parameters that will be adjusted
%injection efficiency
% eta0 = 0.1062;
%electron lifetime, nsec
% taun0 = 0.0355;
%how much faster peak of gain spectra shifts than laser modes 
% Z0 = 5;
%energy of gain spectra peak at J0 injection current density, eV
% E0J0 = 1.1894;

%serialize parameters for use in optimizing function
% x0 = [eta0, taun0, Z0, E0J0];

d.sig = 20e-3;
x0  = [d.eta d.taun d.Z d.E0J0];
x0 = [0.1705    0.058    8    1.18]

%save spectra to files at each injection current
cfg.writeSpectra = 0;

%integration parameters
%range of QD energies in standard deviation units (sigma)
%   1 accounts for 68.7%  of all QDs
%   2 - 95.45%
%   3 - 99.73%
%however, consider that the "action" happens on the very edge 
%   of energy distribution!
cfg.rangeE = 6;
%number of distretization points
cfg.ptsE = 300;


options = optimset('Display', 'iter', 'MaxIter', 50);
exitflag = 0;
while exitflag ~= 1
    x = x0;
    [out, d] = calcX(x, d, cfg);
    [s, powerFine] = powerScale(out, d);
    out.S1Arr = out.S1Arr * s;
    out.S2Arr = out.S2Arr * s;
    disp(x)
    draw(out, d)
    drawnow   
    [x ,fval, exitflag, output] = ...
        fminsearch(@(x) f(x, d, cfg), x0, options);
    x0 = x;
end


function [out, d] = calcX(x, d, cfg)
d.eta = x(1);
d.taun = x(2);
d.Z = x(3);
d.E0J0 = x(4);

JArr = unique(sort([d.exp.JS1x; d.exp.JS2x]));
out = calc(JArr, d, cfg);
end


%penalty function, shows how far model results are from the experiment
%mode intensities are the point of comparison
%the plan is to normalize experiment's and model's mode intensities and
%   interpolate them to a common grid
%the fine at every point is a relative error
function [fine] = f(x, d, cfg)

[out, d] = calcX(x, d, cfg);
[s, fine] = powerScale(out, d);
end