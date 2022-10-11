%{
Task 5 B: Single-Phase Phasor Diagram for the Series RLC Circuit Based on
the measurement in Task 1D, Step 4, draw (using MATLAB or Excel, not by
hand!) the phasor diagram for this case showing the applied voltage V1,
input current I1, inductor voltage V2 and capacitor voltage V3. Then,
repeat the same for the measurement of Step 5. Based on your diagrams, how
does the voltages across the inductor and capacitor relate to the applied
voltage?

%}
% Step 4 is all capacitors, 1 resistor
% Step 5 is all capacitors and all resistors
clear; clc;

V1 = 15*exp(deg2rad(-0.1)*1i);
pfa1 = deg2rad(2.2);
pfa2 = deg2rad(87.2);
pfa3 = deg2rad(-87.5);
I1 = 0.59 * exp((angle(V1)-pfa1)*1i);
V2 = 14.4 * exp((angle(V1)-pfa2)*1i);
V3 = 15.05 * exp((angle(V1)-pfa3)*1i);

% Drawing Phasors
quiver(0,0,real(V1),imag(V1),'Color','b');
xlim([-20 20])
axis equal
grid on;
hold on;
quiver(0,0,real(V2),imag(V2),'Color','r');
quiver(0,0,real(V3),imag(V3),'Color','g');
title('Task 1D Step 4 phasors')
xlabel('real axis (V)');
ylabel('imaginary axis (VR)');
text(real(V1)/2, imag(V1)+2, sprintf('V1(%.2f,%.2f)',[real(V1) imag(V1)]));
text(real(V2), imag(V2), sprintf('V2(%.2f,%.2f)',[real(V2) imag(V2)]));
text(real(V3), imag(V3), sprintf('V3(%.2f,%.2f)',[real(V3) imag(V3)]));
quiver(0,0,real(I1),imag(I1));