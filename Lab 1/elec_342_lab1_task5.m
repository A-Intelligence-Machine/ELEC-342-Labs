%{
Matlab file for ELEC 342 Lab 1 AC/DC Circuits and Basic Measurements
Author: Braedan Chan Date: Sept, 26, 2022
%}
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



%{
%% Task 5.B: Single-Phase Phasor Diagram for the Series RLC circuit
Based on the measurement in Task 1D, Step 4, draw the phasor diagram for
this case showing applied phase voltage V1, input current I1, inductor
voltage V2 and capacitor voltage (V3)
Repeat this process for the sama measurement in Step 5.
Based on your diagrams, how does the voltage across the inductor and
capacitor relate to the applied voltage?
%}

%{
%% Task 5.C: Three-Phase Phasor Diagrams for Wye-connected RL load
Based on the measurement in task 2A, step 3, draw the phasor diagram for
this case showing applied voltages V1, V2, V3, and input currents I1, I2,
I3.
Repeat this process for Step 4 and Step 5
Based on your diagrams, how does the unbalanced load affect the phase
voltages when the neutral point on the load is not connected to the source
neutral?
What is the effect of connecting the neutral points in step 5?
What would happen in the case of an unbalanced load if the load phases are
connected as Delta as in Task 2B?
Explain based on the measurements in Table 6
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

V1a = 10.01*sqrt(2)*sin(2*pi*60/1000*t+3.5*pi/180); % The assumed measured value of V1 if it were an actual sinusoid.
V2a = 10.31*sqrt(2)*sin(2*pi*60/1000*t-238.6*pi/180);
% Plotting out the figure
figure(1); clf;
plot(t,V1,t,V2,t,V3,t,V1a,t,V2a);
title('Task 2A.1 Voltage vs Time');
xlabel('Time (ms)');
ylabel('Voltage (V)');
legend('V1','V2','V3');

%{
figure(2); clf;
plot(t,I1,t,I2,t,I3);
title('Task 2A.2 Current vs Time');
xlabel('Time (ms)');
ylabel('Current (A)');
legend('I1','I2','I3');
%}

figure(2); clf;
vec = [t,V1];
T = (t(2)-t(1))*1e-3; % Set the sample period by the difference between two adjacent time points in ms
Fs = 1/T;
t_vec = 0:1/Fs:t(1,end); % Set up an even time vector by using the sample period and the last time measurement
L = length(V1);
X = fft(V1,L);
X = X(1:L/2); % Take half the frequency spectrum
mx = abs(X); % Magnitude;
f = (0:L/2-1)*Fs/L; % Create the frequency vector
stem(f,mx,Linewidth=2.5);
xlim([0 120]); % Set where you want to look
xlabel("Frequency (Hz)", 'FontSize', 17);
ylabel('Voltage', 'FontSize', 17);
title('Voltage vs. Frequency of Input Voltage from Task 2A, Step 3', 'FontSize', 17);

figure(3); clf;
% Create a phasor plot for the applied phase voltages and currents
Vi = [10.01*(cos(deg2rad(3.5)+i*sin(deg2rad(3.5)))) 10.31*(cos(deg2rad(-238.6))+i*sin(deg2rad(-238.6))) 9.72*(cos(deg2rad(-116.7))+i*sin(deg2rad(-116.7)))];
plotArrow(abs(Vi(1)), angle(Vi(1)), 'r');
plotArrow(abs(Vi(2)), angle(Vi(2)), 'b');
plotArrow(abs(Vi(3)), angle(Vi(3)), 'g');
% And now to plot the current phasors
title('Voltage Phasors Task 5C')

figure(4); clf;
I1i = 0.28 * exp(deg2rad(-41.1)*j);
I2i = 0.29 * exp(deg2rad(-77.5)*j);
I3i = 0.29 * exp(deg2rad(-162.5)*j);
plotPhasor(abs(I1i), angle(I1i));
plotPhasor(abs(I2i), angle(I2i));
plotPhasor(abs(I3i), angle(I3i));


%Draw an Arrowhead and a small circle to show the x-axis intercept).
% Take as input the magnitude and phase in radians. Also the color you want
% it plotted in
function plotArrow(A,phi,c)
hold on

%Define Arrow.
x=[0 A A-0.2 A A-0.2]';
y=[0 0 0.2 0 -0.2]';
%Rotate Arrow.

x1=x*cos(phi)-y*sin(phi);
y1=x*sin(phi)+y*cos(phi);
%Plot Arrow.
plot(x1,y1,c,'LineWidth',2)

%Define small circle and plot it on x axis.
theta=0:pi/4:2*pi;
xp=A*cos(phi)+0.1*cos(theta);
yp=0.1*sin(theta);
fill(xp,yp,c,'EdgeColor',c);

%Draw dotted circle to show magnitude of error.
theta=0:0.1:2*pi;
plot(A*cos(theta),A*sin(theta),char(c,':'));
grid on;
end

function plotPhasor(A,phi)
hold on;
% Starting position of vector
x_start_pos = 0;
y_start_pos = 0;

% Head position of vector
x_final_pos = A*cos(phi);
y_final_pos = A*sin(phi);
% Creating arrow 
hArr = annotation('arrow',[x_start_pos x_final_pos],[y_start_pos y_final_pos]);

end
%{
adft = fft(V2);
V2m = abs(adft(61))
V2p = angle(adft(61))

adft = fft(V3);
V3m = abs(adft(61))
V3p = angle(adft(61))

adft = fft(I1);
I1m = abs(adft(61))
I1p = angle(adft(61))

adft = fft(I2);
I2m = abs(adft(61))
I2p = angle(adft(61))

adft = fft(I3);
I3m = abs(adft(61))
I3p = angle(adft(61))
%}

%{
Art of the fft
Use Fourier transforms to find the frequency components of a signal buried in noise.

Eg. Specify the parameters of a signal with a sampling frequency of 1 kHz and a signal duration of 1.5 seconds.

Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector

Form a signal containing a 50 Hz sinusoid of amplitude 0.7 and a 120 Hz sinusoid of amplitude 1.
S = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t);

Corrupt the signal with zero-mean white noise with a variance of 4.
X = S + 2*randn(size(t));

Plot the noisy signal in the time domain. It is difficult to identify the frequency components by looking at the signal X(t).

Compute the Fourier transform of the signal.
Y = fft(X);

Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

Define the frequency domain f and plot the single-sided amplitude spectrum P1. The amplitudes are not exactly at 0.7 and 1, as expected, 
because of the added noise. On average, longer signals produce better frequency approximations.
f = Fs*(0:(L/2))/L;
plot(f,P1) 
title("Single-Sided Amplitude Spectrum of X(t)")
xlabel("f (Hz)")
ylabel("|P1(f)|")

Now, take the Fourier transform of the original, uncorrupted signal and retrieve the exact amplitudes, 0.7 and 1.0.
Y = fft(S);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

plot(f,P1) 
title("Single-Sided Amplitude Spectrum of S(t)")
xlabel("f (Hz)")
ylabel("|P1(f)|")

%}

% Now the matrix A contains all the information from the file.
% 

%{
%% Task 5.D: Harmonics in AC mains
Use any of the measured/recorded voltages in task 1 or 2, first plot the
recorded voltage waveform.
How close is this waveform to an ideal sinusoid?
Use the MATLAB FFT command to calculate the harmonic constant.
Also plot this content where the harmonic amplitudes are either in volts or
normalized with respect to the fundamental component.
The frequency axis should have units of Hz.
Read the reference material on fft()
You will need to rescale the output of fft to display the harmonic spectrum
properly.
What do you think can cause harmonics in AC systems and in the Lab in
particular?
%}
%{
%% Task 5.E: Single-Phase AC-DC rectifiers
Use the measured/recorded current from task 3A, first plot the recorded
current waveform. How close is this waveform to an ideal sinusoid?
Use fft() to calculate the harmonic content.
Also plot this content where the harmonic amplitudes are either in Amps or
normalized with respect to the fundamental component same as in task 5D
What could be the result of such loads on AC mains?
%}
%{
%% Task 5.F: Three-Phase AC-DC Rectifiers
Use the measured/recorded currents from Task 4, Step 4
Plot the recorded currents I1, I2, I3 waveforms.
Explain the difference between I2 and I3.

Use fft() to calculate the harmonic content of the phase current I1.
How does the waveform of the phase current I1 in Task 4 differ from the
current I1 in Task 3A?
%}