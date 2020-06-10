%draw results of calc()
%   out - calculation results
%   d - laser data
function draw(out, d)

subplot(3,2,1)
plot(out.JArr, out.efArr);
xlabel('J, A/cm^2')
ylabel('eV')
title('Fermi energy')


subplot(3,2,2)
plot(out.JArr, out.nArr, '-o')
ylabel('El. density, 1/cm^3')
yyaxis right
plot(out.JArr, out.S1Arr, 'x-', out.JArr, out.S2Arr, 'o-')
ylabel('Photon density, 1/cm^3')
xlabel('A/cm^2')
legend('n', 'S1', 'S2')


subplot(3,2,3)
plot(out.JArr, out.e1Arr, 'b-x', out.JArr, out.e2Arr, 'r-x', ...
    out.JArr, out.e0Arr, 'g-x')
ylabel('eV')
xlabel('A/cm^2')
legend('Mode 1', 'Mode 2',  'Gain peak')
title('Trans. energy')


subplot(3,2,4)
plot(out.JArr, out.Gain1 - d.alpha1, '-.', out.JArr, out.Gain2 - d.alpha2, '-.')
yyaxis right
plot(out.JArr, out.Rsp1Arr, '-x', out.JArr, out.Rsp2Arr, '-o');
xlabel('A/cm^2')
ylabel('1/ns')
legend('G_1-\alpha_1', 'G_2 - \alpha_2', 'R_{sp1}', 'R_{sp2}')
title('Gains')


subplot(3,2,5)
semilogy(out.JArr, out.S1Arr, 'b-x', out.JArr, out.S2Arr, 'r-x',...
    d.exp.JS1x, d.exp.JS1y, 'bo', d.exp.JS2x, d.exp.JS2y, 'ro')
legend('S_1, model', 'S_2, model', 'S_1, exper.', 'S_2, exper.')
xlabel('J, A/cm^2')
ylabel('Intensity, a.u.')


subplot(3,2,6)
plot(out.JArr, out.S1Arr, 'b-x', out.JArr, out.S2Arr, 'r-x', ...
    d.exp.JS1x, d.exp.JS1y, 'bo', d.exp.JS2x, d.exp.JS2y, 'ro')

legend('S_1, model', 'S_2, model', 'S_1, exper.', 'S_2, exper.')
xlabel('J, A/cm^2')
ylabel('Intensity, a.u.')

end