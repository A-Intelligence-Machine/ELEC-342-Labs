% Task 5: Calculation and Analysis 
clc; clear;

%{
%% Task 5.A: Instantaneous, Real, and Reactive Power in a Parallel RCL Circuit
Calculate the instantaneous, real, and reactive power based on your
recorded data from Task 1C, Step 3
Plot on separate sub-plots so you can clearly see all calculated powers one
after another.
Repeat this process for Task 1C steps 4 and 5
How has adding an inductor and capacitor changed the reactive and real
power?
2. Based on the measurements in Task 1C Step 5, draw the phasor diagram of
this case showing applied voltage V1, Input current I1, Inductor current I2
and Capacitor Current I3
%}

% Given V1 = 20.28 * exp(deg2rad(-0.6)*j)
% Instantaneous Power
%{
P = Vm*sin(w*t+thetaV)*Im*sin(w*t+thetaI)
P = Vm*Im*(sin(w*t+theta_V)*sin(w*t+theta_I)
%}
% Voltage Phasor
V1 = 20.28 * exp(deg2rad(-0.6)*j)
V1p = angle(V1);
V1m = abs(V1);
I1 = 1.02; % Amps rms
pfa = deg2rad(0.7); % IN degrees theta_V-theta_I;
P = 20.6; % In watts the average power
w = 60*2*pi;
t = linspace(0,30e-3);
Pinst = 2 * V1m * I1.*(cos(w.*t+V1p).*cos(w.*t+V1p-pfa));

Preal = V1m * I1*(cos(w*t+V1p).*cos(w*t+V1p-pfa));
Pimag = V1m * I1*sin(w*t+V1p).*sin(w*t+V1p-pfa);

%Pim = imag(V1*conj(I1(exp(j*(V1p-pfa))));

figure(1); clf;
plot(t,Pinst,'color','b');
hold on;
plot(t, Preal, 'color','r');
plot(t,Pimag, 'color','g');
title('Task 5A powers in task 1C One Resistor');
legend('Pinst', 'Preal', 'Pimag');

% R and L
I1 = 1.56;
pfa1 = deg2rad(-45.7);

Pinst = 2 * V1m * I1*(cos(w.*t+V1p).*cos(w.*t+V1p-pfa1));
Preal = V1m * I1*(cos(w*t+V1p).*cos(w*t+V1p-pfa1));
Pimag = V1m * I1*sin(w*t+V1p).*sin(w*t+V1p-pfa1);

%Pim = imag(V1*conj(I1(exp(j*(V1p-pfa))));

figure(2); clf;
plot(t,Pinst,'color','b');
hold on;
plot(t, Preal, 'color','r');
plot(t,Pimag, 'color','g');
title('Task 5A powers in task 1C Resistor+Inductor');
legend('Pinst', 'Preal', 'Pimag')

% R and L and C
I1 = 1.19;
pfa1 = deg2rad(-9.9);
Pinst = 2 * V1m * I1*(cos(w.*t+V1p).*cos(w.*t+V1p-pfa1));
Preal = V1m * I1*(cos(w*t+V1p).*cos(w*t+V1p-pfa1));
Pimag = V1m * I1*sin(w*t+V1p).*sin(w*t+V1p-pfa1);
figure(3); clf;
plot(t,Pinst,'color','b');
hold on;
plot(t, Preal, 'color','r');
plot(t,Pimag, 'color','g');
title('Task 5A powers in task 1C Resistor+Inductor+Capacitor');
legend('Pinst', 'Preal', 'Pimag')

% Drawing out phasors
figure(4);clf;
I1 = 1.19;
pfa1 = deg2rad(-9.9);
I2 = 0.45;
pfa2 = deg2rad(-64.1);
I3 = 1.0;
pfa3 = deg2rad(89.2);
X  = [I1*cos(V1p-pfa1) I2*cos(V1p-pfa2) I3*cos(V1p-pfa3)];
Y = [I1*sin(V1p-pfa1) I2*sin(V1p-pfa2) I3*sin(V1p-pfa3)];
orig = zeros(1,length(X));
q = quiver(orig,orig,X,Y);
axis ([-2 2 -2 2])
grid on;
title('Task 1C RLC current phasors')

text(X(1), Y(1), sprintf('I1(%.2f,%.2f)',[X(1) Y(1)]));
text(X(2), Y(2), sprintf('I2(%.2f,%.2f)',[X(2) Y(2)]));
text(X(3), Y(3), sprintf('I3(%.2f,%.2f)',[X(3) Y(3)]));

figure(5); clf;
quiver(0,0,V1m*cos(V1p), V1m*sin(V1p));
axis([-25 25 -25 25]);
grid on;
text(V1m*cos(V1p)-5, V1m*sin(V1p)-5, sprintf('V1 (%.2f,%.2f)', [V1m*cos(V1p) V1m*sin(V1p)]))
title('Task 1C RLC Voltage Phasor');
