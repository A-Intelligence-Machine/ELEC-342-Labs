% Task 5E
%{
Use the measured/recorded current from Task 3A, first plot the recorded
current waveform. How close is this waveform to an ideal sinusoid? Use
MATLAB FFT commend (or an equivalent tool) to calculate the harmonic
content. Also plot this content where the harmonic amplitudes are either in
Amps or normalized with respect to the fundamental component, same as you
did in Task 5D.

What could be the result of such loads on the AC Mains?

%}
fileID = fopen('lab1D_task3_heavyLoad_BandP.txt','r');
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

I1a = sqrt(2)*1.19*cos(2*pi*60/1000*t+deg2rad(-94.3));

figure(1); clf;
plot(t,I1,t,I1a);
title('Task 3A Current vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');
legend('I1','I1calc');
grid on;

figure(2); clf;
vec = [t,I1];
T = (t(2)-t(1))*1e-3; % Set the sample period by the difference between two adjacent time points in ms
Fs = 1/T;
t_vec = 0:1/Fs:t(1,end); % Set up an even time vector by using the sample period and the last time measurement
L = length(I1);
X = fft(I1,L);
X = X(1:L/2); % Take half the frequency spectrum
mx = abs(X/L); % Magnitude;
f = (0:L/2-1)*Fs/L; % Create the frequency vector
stem(f,2*mx,Linewidth=2.5);
xlim([0 500]); % Set where you want to look
xlabel("Frequency (Hz)");
ylabel('Current');
title('Current vs. Frequency Task 3A');
