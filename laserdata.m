function d = laserdata()
%physical constants

%Planck constant, eV*sec -> eV*nsec
h = 4.135e-15 * 1e9; 
%electron charge, coulomb
d.e = 1.6e-19; 
%speed of light in vacuum, cm/sec -> cm/nsec
c = 29979245800 / 1e9; 


%laser parameters
d.D = 10 * 1e-4; %diameter, mum -> cm
d.area = pi*(d.D/2)^2;

%injection efficiency
% d.eta = 10e-2;
d.eta = 0.0908; %automatic adjustment result
%electron lifetime, nsec
%d.taun = 0.037;
d.taun = 0.0826; %automatic adjustment result
%active area effective thickness, nm -> cm
d.d = 50 / 1e7;
%QD density, 1/cm^2 -> 1/cm^3
d.nQD = 5e11 / d.d;
%how much faster peak of gain spectra shifts than laser modes 
%d.Z = 5.3;
d.Z = 6.3545; %from adjust()
%inhomogeneous broadening, eV
d.sig = 25e-3;
%homogeneous broadening, eV
d.gamma = 5e-3;
%temperature, eV
d.kT = 23e-3;
%maximum gain, 1/nsec
d.Gmax = 510;
%refractive index (assuming GaAs)
d.n = sqrt(12.9);
%light attenuation, 1/nsec
d.alpha = 45;
%inj. current density at reference point, A/cm^2
d.J0 = 4000;
%gain peak energy at reference point, eV
% d.E0J0 = 1.96;
d.E0J0 = 1.1905; %automatic adjustment result

%mode volume, cm^3
d.Vmode = d.area*d.d*d.eta;
%attenuation of both modes, 1/nsec
d.alpha1 = d.alpha;
d.alpha2 = d.alpha;


%1st mode experiment data
%injection current, mA -> A
I1 = [2.8 3 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4 4.2 4.4 4.6 4.8 5].' / 1e3;
%injection current density, A/cm^2
J1 = I1 / d.area;
d.exp.JEx = J1;
%1st mode's energy, 1/cm -> eV
E1 = [9123.7 9122.6 9121.4 9120.8 9120.2 9119.7 9118.9 9118.3 9117.7 ...
    9117 9116.3 9114.67 9113.29 9111.81 9110.39 9108.89].' * h*c;
d.exp.JEy = E1;
d.exp.JS1x = J1;
%1st mode intensity
d.exp.JS1y = [0.0098 0.0163 0.0249 0.0336 0.0482 0.0917 0.158 0.31 0.631...
    0.789 0.96978 0.52479 0.1451 0.08283 0.05986 0.04776].';
%interpolating to common current density array
d.E1J0 = interp1(J1, E1, d.J0, 'linear', 'extrap');

p = polyfit(J1, E1, 1);
d.chi = p(1); %eV/(A/cm^2)


%2nd mode experiment data
%injection current, mA -> A
I2 = [3.4 3.6 3.8 4 4.1 4.2 4.4 4.6 4.8 4.9 5 5.2].' / 1e3;
%injection current density, A/cm^2
J2 = I2 / d.area;
d.exp.JE2x = J2;
%1/cm -> eV
%mode energy, 1/cm -> eV
E2 = h*c*[9032.416 9031.158 9029.908 9028.63 9027.995 9027.307 9026.4 ...
    9025.03 9023.65 9022.76 9022.14 9020.27].';
d.exp.JE2y = E2;
%interpolating to common current density array
d.E2J0 = interp1(J2, E2, d.J0, 'linear', 'extrap');

%injection current, mA -> A
I2 = [3.4 3.6 3.8 4 4.1 4.2 4.4 4.6 4.8 4.9 5 5.2 3.3 3.2 3 2.8].' / 1e3;
%injection current density, A/cm^2
J2 = I2 / d.area;
d.exp.JS2x = J2;
%mode intensity
S2 = [0.0315 0.0648 0.137 0.259 0.403 1.54 4.67 7.05 9.26 9.58 10.05 ...
    10.8 0.0235 0.0164 0.0102 0.0062].';
d.exp.JS2y = S2;

end