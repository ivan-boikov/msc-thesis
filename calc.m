
function out = calc(JArr, laser, cfg)

N = length(JArr);

%preallocating output arrays
zeroN = zeros(N, 1);
%current density, A/cm^2
out.JArr = JArr;
%bound electron density, 1/cm^3
out.nArr = zeroN;
%Fermi energy, eV
out.efArr = zeroN;
%photon density in modes 1 and 2, 1/cm^3
out.S1Arr = zeroN;
out.S2Arr = zeroN;
%modal gain of modes 1 and 2, 1/nsec
out.Gain1 = zeroN;
out.Gain2 = zeroN;
%spontaneous emission rate in modes 1 and 2, 1/nsec
out.Rsp1Arr = zeroN;
out.Rsp2Arr = zeroN;
%energies of gain peak and modes 1 and 2, eV
out.e0Arr = zeroN;
out.e1Arr = zeroN;
out.e2Arr = zeroN;
%energy of gain spectra peak, eV
out.GmaxE = zeroN;
%gain spectra FWHM, eV
out.GainFWHM = zeroN;

for i = 1:N
    J = JArr(i);
    
    %starting point, does not have to be close to solution, 
    %   just being sensible is enough
    ef = laser.E0J0;
    S1 = laser.nQD;
    S2 = laser.nQD;
    
    %integrating the system until equilibrium
    x0 = [ef; S1; S2];
    [~, x] = ode15s(@(t,x) deriv(x, J, laser, cfg), [0, cfg.integT], x0);
    x = x(end, :);
    ef = x(1);
    S1 = x(2);
    S2 = x(3);
        
    r = interParam(x, J, laser, cfg);
    
    if cfg.writeSpectra
        %used for drawing plots in Origin
        csvwrite(sprintf('data/spectra%.0fk.csv', JArr(i) / 1e3), ...
            [r.E, r.rho/max(r.rho), r.f, r.F1, r.F2, r.G/laser.Gmax])
        
        %finding energy of gain spectra peak
        [~,ii] = max(r.G);
        %max() is not enough - there will be "stairs" due to dicretization
        %   interpolating for better precision
        F = griddedInterpolant(r.E, r.G, 'spline');
        out.GmaxE(i) = arrayfun(@(xx)fminsearch(@(x)-F(x), xx), r.E(ii));
        
        out.GainFWHM(i) = fwhm(r.E, r.G);
    end
    
    out.efArr(i) = ef;
    out.S1Arr(i) = S1;
    out.S2Arr(i) = S2;
    out.nArr(i) = r.nEl;
    out.Gain1(i) = r.G1;
    out.Gain2(i) = r.G2;
    out.Rsp1Arr(i) = r.Rsp1;
    out.Rsp2Arr(i) = r.Rsp2;
    out.e0Arr(i) = r.e0;
    out.e1Arr(i) = r.e1;
    out.e2Arr(i) = r.e2;    
end

end


%evaulate derivatives of Fermi energy and photon densities
%   x - serialized parameters
%   J - injection current density, coulomb/nsec/cm^2
function der = deriv(x, J, laser, cfg)
S1 = x(2);
S2 = x(3);

par = interParam(x, J, laser, cfg);

der = x*0;

dndEf = 2*laser.nQD * trapz(par.rho.*par.f) * par.dE;

der(1) = laser.eta*(J/1e9)/(laser.e*laser.d) - par.nEl/laser.taun - ...
    - par.Rsp1/laser.Vmode - par.G1*S1 - par.Rsp2/laser.Vmode - par.G2*S2;
der(1) = der(1) / dndEf;
der(2) = (par.G1 - laser.alpha1)*S1 + par.Rsp1/laser.Vmode;
der(3) = (par.G2 - laser.alpha2)*S2 + par.Rsp2/laser.Vmode;

end



%calculate intermediate laser parameters
%   x - serialized Fermi energy and photon densities
%   J - injection current density, coulomb/nsec/cm^2
%   d - laser data
%   cfg - configuration parameters
function out = interParam(x, J, d, cfg)

Ef = x(1);

dJ = J-d.J0;
%energies of gain spectra peak and two modes
e0 = d.E0J0 + d.chi*d.Z*dJ;

E = linspace(e0-d.sig*cfg.rangeE, e0+d.sig*cfg.rangeE, cfg.ptsE).';
dE = E(2) - E(1);

%QD energy distribution
rho = 1./sqrt(2*pi)/d.sig * exp(-(E - e0).^2 / (2*d.sig^2));
%fill-factor
f = 1./(1 + exp((E - Ef)/d.kT));

%bound electron density
out.nEl = 2*d.nQD * trapz(rho.*f) * dE;


%mode 1
%energy
e = d.E1J0 + d.chi*dJ;
out.e1 = e;
%homogeneous broadening
F = 1./(1 + (E - e).^2/d.gamma^2);
denom = trapz(rho.*F);
%modal gain
out.G1 = d.Gmax * trapz(rho.*F.*(2*f-1)) / denom;
%sponateous emission rate
out.Rsp1 = d.Gmax * trapz(rho.*F.*f) / denom;


%same, mode 2
e = d.E2J0 + d.chi*dJ;
out.e2 = e;
F = 1./(1 + (E - e).^2/d.gamma^2);
denom = trapz(rho.*F);
out.G2 = d.Gmax * trapz(rho.*F.*(2*f-1)) / denom;
out.Rsp2 = d.Gmax * trapz(rho.*F.*f) / denom;

%these are needed regardless you want to write spectra or not
out.e0 = e0;
out.E = E;
out.dE = dE;
out.rho = rho;
out.f = f;

if cfg.writeSpectra
    out.F1 = F;
    out.F2 = F;
    out.G = d.Gmax * (rho.*(2*f-1)) / max(rho);
end

end



%calculate full width at half maximum assuming peak is single
%   x - coordinates
%   y - values
%   w - FWHM
function width = fwhm(x, y)

[ymax, imax] = max(y);

il = imax;
while y(il) > ymax/2
    il = il - 1;
end
xleft = interp1([y(il) y(il+1)], [x(il) x(il+1)], ymax/2);

ir = imax;
while y(ir) > ymax/2
    ir = ir + 1;
end
xright = interp1([y(ir-1) y(ir)], [x(ir-1) x(ir)], ymax/2);

width = xright - xleft;

end