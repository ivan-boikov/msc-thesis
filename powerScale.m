%get scaling factor of photon densities to compare them with the
%   experiment
%since the output power from experiment is in a.u., we can multiply 
%   photon density by any number
%   out - output structure of calc() function
%   d - laser data structure from laserdata() function
%   s - scaling factor for photon densities
%   fine - how far scaled photon densities deviate from the experiment
function [s, fine] = powerScale(out, d)  
    %
    f = @(s) arrFine(s, out.JArr, out.S1Arr, d.exp.JS1x, d.exp.JS1y) +...
        5*arrFine(s, out.JArr, out.S2Arr, d.exp.JS2x, d.exp.JS2y);
    
     options = optimset('Display', 'none', 'TolX', 1e-8, 'MaxIter', 5000);
     
    [s, fine, exitflag] = fminsearch(f, 1, options);
    
    %if argument/function tolerance was met
    if exitflag ~= 1
        disp('WARNING: scale() failed to converge')
    end
end

%calculate how far (x,y) curve is from (xTgt, yTgt) based on relative error
function res = arrFine(s, x, y, xTgt, yTgt)
    %interpolate to a common grid
    yy = interp1(x, y, xTgt) * s;
    %average sum of relative errors per point - if 2nd mode has more ...
    %data points, we don't want it to have an advantage over 1st mode
    res = norm((yTgt - yy)./yTgt) / length(xTgt);
end