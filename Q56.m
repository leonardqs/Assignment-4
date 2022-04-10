clear
close all

R1 = 1;
Cap = 0.25;
R2 = 2;
Ind = 0.2;
R3 = 20; % Assume R3 = 20 in Part 5 and 6
alpha = 100;
R4 = 0.1;
Ro = 1000;
Cn = 0.00001; %0.00001 %0.001 %0.1 %10
% Vout with difeerent Cn
G1 = 1/R1;
G2 = 1/R2;
G3 = 1/R3;
G4 = 1/R4;
Go = 1/Ro;

G = zeros(7,7);
C = zeros(7,7);

G(1,1) = G1;
G(1,2) = -G2-G1;
G(1,3) = -G3;
G(2,2) = 1;
G(2,3) = -1;
G(3,3) = -G3;
G(3,6) = 1;
G(4,4) = 1;
G(4,6) = -alpha;
G(5,4) = G4;
G(5,5) = -G4-Go;
G(6,1) = 1;
G(7,1) = -G1;
G(7,2) = G1;
G(7,7) = 1;


C(1,1) = Cap;
C(1,2) = -Cap;
C(2,6) = -Ind;
C(3,3) = -Cn;
C(7,1) = -Cap;
C(7,2) = Cap;


TimeStep = 1e-3; % 1e-3 % 1e-4 %1e-2
% Vout with different time steps 

In = 0.001*randn(1/TimeStep,1); % In=0.01 to make noise visible


time = linspace(0,1,1/TimeStep);
figure(1)
subplot(311)
Input1 =  ones(length(time),1);
Input1(1,1) = 0;
Input1(2,1) = 0;
Input1(3,1) = 0;
Input1(4,1) = 0;
plot(time,Input1)
title('Input Signal 1');
xlabel('Time (s)')
ylabel('Input signal amplitude (V)')

subplot(312)
f2 = 1/0.03; % Input Signal 2 Frequency
Input2 = sin(2*pi*f2*time);
plot(time,Input2)
title('Input Signal 2');
xlabel('Time (s)')
ylabel('Input signal amplitude (V)')

subplot(313)

Input3 = exp(-0.5*((time - 0.06).^2) /(0.03)^2);
plot(time,Input3)
xlim([0 0.2])
title('Input Signal 3');
xlabel('Time (s)')
ylabel('Input signal amplitude (V)')

A = C/TimeStep + G;
InvA = inv(A);

V1 = zeros(7,1);
V2 = zeros(7,1);
V3 = zeros(7,1);
F1 = zeros(7,1);
F2 = zeros(7,1);
F3 = zeros(7,1);


Vout1 = zeros(1/TimeStep,1);
Vout2 = zeros(1/TimeStep,1);
Vout3 = zeros(1/TimeStep,1);
for i = 1:length(time)
    
    LastV1 = V1;
    LastV2 = V2;
    LastV3 = V3;
    
    F1(6) = Input1(i);
    F2(6) = Input2(i);
    F3(6) = Input3(i);

    F1(3) = In(i);
    F2(3) = In(i);
    F3(3) = In(i);
    

    V1 = A\(C*(LastV1/TimeStep) + F1);
    Vout1(i,1) = V1(5,:);
    V2 = A\(C*(LastV2/TimeStep) + F2);
    Vout2(i,1) = V2(5,:);
    V3 = A\(C*(LastV3/TimeStep) + F3);
    Vout3(i,1) = V3(5,:);
end

figure(2)
subplot(311)
plot(time,Vout1)
title('Output Signal 1');
xlabel('Time (s)')
ylabel('Output signal amplitude (V)')
subplot(312)
plot(time,Vout2)
title('Output Signal 2');
xlabel('Time (s)')
ylabel('Output signal amplitude (V)')
subplot(313)
plot(time,Vout3)
title('Output Signal 3');
xlabel('Time (s)')
ylabel('Output signal amplitude (V)')



% Sample frequency = 1000/s

frequency = linspace(-500,500,1/TimeStep); % FFT Frequency Vector 

ftIn1 = fft(Input1);
shiftIn1 = fftshift(ftIn1);
ftIn2 = fft(Input2);
shiftIn2 = fftshift(ftIn2);
ftIn3 = fft(Input3);
shiftIn3 = fftshift(ftIn3);

ftOut1 = fft(Vout1);
shiftOut1 = fftshift(ftOut1);
ftOut2 = fft(Vout2);
shiftOut2 = fftshift(ftOut2);
ftOut3 = fft(Vout3);
shiftOut3 = fftshift(ftOut3);

figure(3)
subplot(311)
plot(frequency,abs(ftIn1));
title('FFT Input Signal 1');
xlabel('Frequency (rad/s)')
ylabel('Input signal amplitude (V)')
subplot(312)
plot(frequency,abs(ftIn2));
title('FFT Input Signal 2');
xlabel('Frequency (rad/s)')
ylabel('Input signal amplitude (V)')
subplot(313)
plot(frequency,abs(ftIn3));
title('FFT Input Signal 3');
xlabel('Frequency (rad/s)')
ylabel('Input signal amplitude (V)')

figure(4)
subplot(311)
plot(frequency,abs(shiftIn1));
title('FFT Shift Input Signal 1');
xlabel('Frequency (rad/s)')
ylabel('Input signal amplitude (V)')
subplot(312)
plot(frequency,abs(shiftIn2));
title('FFT Shift Input Signal 2');
xlabel('Frequency (rad/s)')
ylabel('Input signal amplitude (V)')
subplot(313)
plot(frequency,abs(shiftIn3));
title('FFT Shift Input Signal 3');
xlabel('Frequency (rad/s)')
ylabel('Input signal amplitude (V)')

figure(5)
subplot(311)
plot(frequency,abs(ftOut1));
title('FFT Output Signal 1');
xlabel('Frequency (rad/s)')
ylabel('Output signal amplitude (V)')
subplot(312)
plot(frequency,abs(ftOut2));
title('FFT Output Signal 2');
xlabel('Frequency (rad/s)')
ylabel('Output signal amplitude (V)')
subplot(313)
plot(frequency,abs(ftOut3));
title('FFT Output Signal 3');
xlabel('Frequency (rad/s)')
ylabel('Output signal amplitude (V)')

figure(6)
subplot(311)
plot(frequency,abs(shiftOut1));
title('FFT Shift Output Signal 1');
xlabel('Frequency (rad/s)')
ylabel('Output signal amplitude (V)')
subplot(312)
plot(frequency,abs(shiftOut2));
title('FFT Shift Output Signal 2');
xlabel('Frequency (rad/s)')
ylabel('Output signal amplitude (V)')
subplot(313)
plot(frequency,abs(shiftOut3));
title('FFT Shift Output Signal 3');
xlabel('Frequency (rad/s)')
ylabel('Output signal amplitude (V)')



