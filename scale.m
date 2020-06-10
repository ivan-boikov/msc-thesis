function s = scale(out, d)  
    f = @(s) fine(s, out.JArr, out.S1Arr, d.exp.JS1x, d.exp.JS1y) +...
        5*fine(s, out.JArr, out.S2Arr, d.exp.JS2x, d.exp.JS2y);
    options = optimset('Display', 'iter', 'TolX', 1e-8, 'MaxIter', 5000);
    [s, ~, exitflag] = fminsearch(f, 1, options);
    if exitflag ~= 1
        disp('WARNING: scale() failed to converge')
    end
end

function res = fine(s, x, y, xTgt, yTgt)
    yy = interp1(x, y, xTgt) * s;   
    res = norm((yTgt - yy)./yTgt) / length(xTgt);
end