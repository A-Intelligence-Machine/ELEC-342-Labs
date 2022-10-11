%{
Task 5 F: Three-Phase AC-DC Rectifiers Use the measured/recorded currents
from Task 4, Step4, and plot the recorded currents I1, I2, I3 waveforms.
Explain the difference between I2 and I3.

Use MATLAB FFT command to calculate the harmonic content of the phase
current I1. How does the waveform of the phase current I1 in Task 4 is
different from current I1 in Task 3A

%}

fileID = fopen('lab1D_task4_rect_heavyLoad_BandP.txt','r');
fgetl(fileID); % Read and discard the first line

[A count] = fscanf(fileID,'%f %f %f %f %f %f %f\n',[7 Inf]);

fclose(fileID);
t = A(1,:);
V1 = A(2,:);
V2 = A(3,:);
V3 = A(4,:);
I1 = A(5,:);
I2 = A(6,:);
I3 = A(7,:);

figure(1); clf;
plot(t,I1,t,I2,t,I3);
legend('I1','I2','I3');
xlim([0 30]); % Set where you want to look
title('Task 4 Currents vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');

figure(3); clf;
I1a = sqrt(2)*1.39*cos(2*pi*60/1000*t+deg2rad(-79.1));
plot(t,I1,t,I1a);
grid on;
legend('I1','I1calc');
title('Task 4 Current 1 (and assumed sinusoid) vs Time');
xlim([0 50]);
xlabel('Time (ms)');
ylabel('Current (A)');


figure(2); clf;
vec = [t,I1];
T = (t(2)-t(1))*1e-3; % Set the sample period by the difference between two adjacent time points in ms
Fs = 1/T;
t_vec = 0:1/Fs:t(1,end); % Set up an even time vector by using the sample period and the last time measurement
L = length(I1);
X = fft(I1,L);
X = X(1:fix(L/2)); % Take half the frequency spectrum
mx = abs(X/L); % Magnitude;
f = (0:L/2-1)*Fs/L; % Create the frequency vector
stem(f,2*mx,Linewidth=2.5);
xlim([0 500]); % Set where you want to look
xlabel("Frequency (Hz)");
ylabel('Current (A)');
title('Current vs. Frequency Task 4');
