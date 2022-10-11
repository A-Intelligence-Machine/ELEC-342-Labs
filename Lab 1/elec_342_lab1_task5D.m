% Task 5D 
%{
Task 5 D: Harmonics in AC mains Use any of the measured/recorded voltages
in Task 1 or 2, first plot the recorded voltage waveform. How close is this
waveform to an ideal sinusoid? Use MATLAB FFT command (or an equivalent
tool) to calculate the harmonic content. Also plot this content where the
harmonic amplitudes are either in volts or normalized with respect to the
fundamental component. The frequency axis should have units of Hz. You can
read help on FFT on the MATLAB command line and/or read online
documentation to learn how this function should be properly used. You will
need to rescale the output of the FFT to display the harmonic spectrum
properly.
%}


fileID = fopen('lab1D_task2a_BandP.txt','r');
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

% The assumed measured value of V1 if it were an actual sinusoid.
V1a = 10.01*sqrt(2)*sin(2*pi*60/1000*t+3.5*pi/180); 

figure(1); clf;
plot(t,V1,t,V1a);
title('Task 2A Voltage#1 vs Time');
xlabel('Time (ms)');
ylabel('Voltage (V)');
legend('Vmeasured','Vcalculated');
grid on;

figure(2); clf;
vec = [t,V1];
T = (t(2)-t(1))*1e-3; % Set the sample period by the difference between two adjacent time points in ms
Fs = 1/T;
t_vec = 0:1/Fs:t(1,end); % Set up an even time vector by using the sample period and the last time measurement
L = length(V1);
X = fft(V1,L);
X = X(1:L/2); % Take half the frequency spectrum
mx = abs(X/L); % Magnitude;
f = (0:L/2-1)*Fs/L; % Create the frequency vector
stem(f,mx,Linewidth=2.5);
xlim([0 500]); % Set where you want to look
xlabel("Frequency (Hz)");
ylabel('Voltage');
title('Voltage vs. Frequency of Input Voltage from Task 2A, Step 3');