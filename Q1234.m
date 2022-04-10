
clear
close all


VoL = 0.1:0.1:10; % Left side of the Area has Boundary Voltage = VoL


VoR = 0; % Right side of the Area has Boundary Voltage = VoR
VoT = 0; % Top side of the Area has Boundary Voltage = VoT
VoB = 0; % Bottom side of the Area has Boundary Voltage = VoB
Pixel = 5e8; % Number of mesh per unit length/width



Length = 2e-7; % Area Length
Width = 1e-7; % Area Width
Lb = 0.4e-7; % Box Length
Wb = 0.3e-7; % Box Width
xBox = Lb*Pixel; % Box length Pixel Number
yBox = Wb*Pixel; % Box width Pixel Number

BottleNeck =  Width - Wb*2;

nx = Length*Pixel; % Area length Pixel Number
ny = Width*Pixel; % Area wdith Pixel Number
Sigma = 1; % Outside Box Area Conductivity
BoxSigma = 0.01; % Inside Box Area Conductivity
G = sparse(nx*ny,nx*ny); % G matrix has size(nx*ny,nx*ny)
B = zeros(nx*ny,1); % B is the product of G matrix * V
Conductivity = Sigma*ones(nx,ny); % Conductivity of the entire area

q_0 = 1.60217653e-19;             % electron charge
hb = 1.054571596e-34;             % Dirac constant
h = hb * 2 * pi;                % Planck constant
m_0 = 9.10938215e-31;             % electron mass
kb = 1.3806504e-23;               % Boltzmann constant
eps_0 = 8.854187817e-12;          % vacuum permittivity
mu_0 = 1.2566370614e-6;           % vacuum permeability
c = 299792458;                    % speed of light
g = 9.80665;                      % metres (32.1740 ft) per s²



timeStep = 5e-15;                     % Time Step
totalParticle = 1e4;
displayParticle = 10;              % Display Particle in the box
SimulationTime = 300;                       % Simulation Time
T = 300;                            % Default Temperature
t_mn = 0.2e-12;                     % Mean time between collision 
loop_index = 10000;
time = 0;


meanVx = zeros(length(VoL),1);

% exponential scattering probability
Pscat = 1- exp(-(timeStep/t_mn));

for u = 1:1:length(VoL)
    100*u/length(VoL)
    for iRow = 1:nx
        for jColumn = 1:ny
    
            if iRow>=(nx-xBox)/2 && iRow<=((nx+xBox)/2) && jColumn<=yBox
                Conductivity(iRow,jColumn) = BoxSigma;
            elseif iRow>=(nx-xBox)/2 && iRow<=((nx+xBox)/2) && jColumn<=ny && jColumn>ny-yBox
                Conductivity(iRow,jColumn) = BoxSigma;
            end
        end
    end



    for iRow = 1:nx
        for jColumn = 1:ny
            n = jColumn+(iRow-1)*ny;
            % Left side Boundary Condition
           
            if iRow == 1      
                G(n,:) = 0;
                G(n,n) = 1;
                B(n) = VoL(u);
            % Right side Boundary Condition
               
            elseif iRow == nx
                G(n,:) = 0;
                G(n,n) = 1;
                B(n) = VoR;
    
            % Bottom side Boundary Condition
            elseif jColumn == 1    
                nxm = jColumn+((iRow-1)-1)*ny;
                nxp = jColumn+((iRow+1)-1)*ny;
                nyp = (jColumn+1)+(iRow-1)*ny;
    
                rxm = (Conductivity(iRow,jColumn) + Conductivity(iRow-1,jColumn))/2;
                rxp = (Conductivity(iRow,jColumn) + Conductivity(iRow+1,jColumn))/2;
                ryp = (Conductivity(iRow,jColumn) + Conductivity(iRow,jColumn+1))/2;
                
                G(n,n) = -(rxm+rxp+ryp);
                G(n,nxm) = rxm;
                G(n,nxp) = rxp;
                G(n,nyp) = ryp;
                
            % Top side Boundary Condition
            elseif jColumn == ny 
                nxm = jColumn+((iRow-1)-1)*ny;
                nxp = jColumn+((iRow+1)-1)*ny;
                nym = (jColumn-1)+(iRow-1)*ny;
    
                rxm = (Conductivity(iRow,jColumn) + Conductivity(iRow-1,jColumn))/2;
                rxp = (Conductivity(iRow,jColumn) + Conductivity(iRow+1,jColumn))/2;
                rym = (Conductivity(iRow,jColumn) + Conductivity(iRow,jColumn-1))/2;
                
                G(n,n) = -(rxm+rxp+rym);
                G(n,nxm) = rxm;
                G(n,nxp) = rxp;
                G(n,nym) = rym;
    
            else 
                nxm = jColumn+((iRow-1)-1)*ny;
                nxp = jColumn+((iRow+1)-1)*ny;
                nym = (jColumn-1)+(iRow-1)*ny;
                nyp = (jColumn+1)+(iRow-1)*ny;
    
                rxm = (Conductivity(iRow,jColumn) + Conductivity(iRow-1,jColumn))/2;
                rxp = (Conductivity(iRow,jColumn) + Conductivity(iRow+1,jColumn))/2;
                rym = (Conductivity(iRow,jColumn) + Conductivity(iRow,jColumn-1))/2;
                ryp = (Conductivity(iRow,jColumn) + Conductivity(iRow,jColumn+1))/2;
                
                G(n,n) = -(rxm+rxp+rym+ryp);
                G(n,nxm) = rxm;
                G(n,nxp) = rxp;
                G(n,nym) = rym;
                G(n,nyp) = ryp;
            end
    
        end
    
    end





    Vn = G\B; % Find (ny*ny:1) size, V = G\B
    % Mapping the V to the matrix size of nx*ny
    for iRow = 1: nx
        for jColumn = 1:ny
             n = jColumn+(iRow-1)*ny;
             V(iRow,jColumn) = Vn(n);    
        end
    end
    
    x = linspace(0,Length,nx);
    y = linspace(0,Width,ny);
    [X,Y] = meshgrid(x,y);




    % 3D Plot of Electric Field E(x,y)
    [Ex,Ey] = gradient(-V');%.*(1/timeStep);
    ExEff = Ex.*Pixel;
    EyEff = Ey.*Pixel;

    cond = Conductivity';
    
    Jx = cond.*(ExEff);
    Jy = cond.*(EyEff);
    
    mass_eff =  0.26*m_0;
    initial_Vth = sqrt((2* kb *T)/mass_eff);
    particle_Time = zeros(1,totalParticle);
    
    
    x_pos = zeros(1,totalParticle);
    y_pos = zeros(1,totalParticle);
    
    
    index = 1;
    for k = 1:loop_index
        rand_x = 2*rand;
        rand_y = rand;
        if(rand_x < 0.8 || rand_x > 1.2) || ((rand_y > 0.4 && rand_y < 0.6) && (rand_x >= 0.8 || rand_x <= 1.2))
            x_pos(1,index) = 1e-7*rand_x;
            y_pos(1,index) = Width*rand_y;
            index = index + 1;
        end
        if index == totalParticle + 1
            break
        end
    end    


    Vx = initial_Vth.*randn(1,totalParticle);
    Vy = initial_Vth.*randn(1,totalParticle);


    xm = linspace(0,Length,100);
    ym = linspace(0,Width,50);
    [XM,YM] = meshgrid(xm,ym);

    for k = 1: SimulationTime
        
        totalTime = timeStep*SimulationTime;
        % Update the electron position and velocity each time step
        prev_time = time;
        time = prev_time + timeStep;
        %prev_tmn = tmn;
        LastX = x_pos;
        LastY = y_pos;
    %    
        Ex_p = interp2(XM,YM,ExEff,x_pos,y_pos);
        Ey_p = interp2(XM,YM,EyEff,x_pos,y_pos);
    
        XAcceleration = (Ex_p*q_0)/m_0;
        YAcceleration = (Ey_p*q_0)/m_0;
    
        Sx = Vx.*timeStep+ (1/2)*timeStep*timeStep*XAcceleration;
        Vx = Vx + timeStep*XAcceleration;
        Sy = Vy.*timeStep + (1/2)*timeStep*timeStep*YAcceleration;
        Vy = Vy + timeStep*YAcceleration;
        
        x_pos = LastX + Sx;
        y_pos = LastY + Sy;
        
        
    
        VTot=sqrt(Vx.^2 + Vy.^2);
        LastX(x_pos<0) = Length;
        x_pos(x_pos<0)= x_pos(x_pos<0)+Length;
    
        
        LastX(x_pos>Length) = 0;
        x_pos(x_pos>Length)=x_pos(x_pos>Length)-Length;
    
        Vy(y_pos<=0) = -Vy(y_pos<=0);
        y_pos(y_pos<=0)=0;
        
        Vy(y_pos>=Width) = -Vy(y_pos>=Width);
        y_pos(y_pos>=Width)=Width;
        shouldScat =  Pscat>rand(n,1);
        particle_Time(Pscat>rand(n,1))=0;
        particle_Time = particle_Time + timeStep;
        Vx(shouldScat)=initial_Vth/sqrt(2)*randn(sum(shouldScat),1);%(1,sum(shouldScat));
        Vy(shouldScat)=initial_Vth/sqrt(2)*randn(sum(shouldScat),1);%(1,sum(shouldScat));
        VTot=sqrt(Vx.^2 + Vy.^2);
        for n = 1:totalParticle
    
          
            if (x_pos(1,n) >= 0.8e-7 && x_pos(1,n) <= 1e-7 ) && (y_pos(1,n) <= 0.4e-7 || y_pos(1,n) >= 0.6e-7) && LastX(1,n) < 0.8e-7
                x_pos(1,n) = 0.8e-7;
                LastX(1,n) = 0.8e-7;
                Vx(1,n) = -(Vx(1,n));
        
            elseif (x_pos(1,n) >= 1e-7 && x_pos(1,n) <= 1.2e-7) && (y_pos(1,n) <= 0.4e-7 || y_pos(1,n) >= 0.6e-7) && LastX(1,n) > 1.2e-7  
                x_pos(1,n) = 1.2e-7;
                LastX(1,n) = 1.2e-7;
                Vx(1,n) = -(Vx(1,n));
        
            elseif (y_pos(1,n) <= 0.4e-7) && (x_pos(1,n) >= 0.8e-7 && x_pos(1,n) <= 1.2e-7)
                y_pos(1,n) = 0.4e-7;
                LastY(1,n) = 0.4e-7;
                Vy(1,n) = -(Vy(1,n));
            elseif (y_pos(1,n) >= 0.6e-7) && (x_pos(1,n) >= 0.8e-7 && x_pos(1,n) <= 1.2e-7)
                y_pos(1,n) = 0.6e-7;
                Vy(1,n) = -(Vy(1,n));
            end
        
        end

    
    end
    meanVx(u) = mean(Vx);

end




figure(1)
Ix = q_0*((1e19)*(2e-7))*meanVx; % Current = Mean Drift Velocity * Charge
plot(VoL, Ix);
title('Current Ix vs Average Voltage');
xlabel('Average Voltage (V)')
ylabel('Current Ix (A)')


bestfit = polyfit(VoL,Ix,1);

Rv = 1/bestfit(1,1);


R1 = 1;
Cap = 0.25;
R2 = 2;
L = 0.2;
R3 = Rv; 
% R3 Value is the linear fit of Current vs Voltage to determine 
% the resistance value of the device and use this value as R3
alpha = 100;
R4 = 0.1;
Ro = 1000;

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
G(3,3) = G3;
G(3,6) = -1;
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
C(2,6) = -L;
C(7,1) = -Cap;
C(7,2) = Cap;



freq = 0;
% DC case sweep the input voltage V1 from -10V to 10V and plot VO
and the voltage at V3.
VinDC = [-10:1:10];
VinDC = VinDC';
VoutDC = zeros(length(VinDC),1);
V3DC = zeros(length(VinDC),1);
F = zeros(7,1);
for i = 1:1:length(VinDC)
    F(6) = VinDC(i);
    V_DC =(G+(1i*2*pi*freq*C))\F; 
    VoutDC(i,1)= V_DC(5,1);  
    V3DC(i,1)= V_DC(3,1); 
end
figure(2)
plot(VinDC,VoutDC)
hold on
plot(VinDC,V3DC)
legend('Vout','V3');
title('DC Voltages vs DC Input Voltage');
xlabel('DC Input Voltage (V)')
ylabel('DC Voltage (V)')

freq = linspace(0,100,100);
VoutAC = zeros(length(freq),1);
for i = 1:1:100
%     Vout = V3
    F(6) = 1; % set input Voltage = 1
    k = (G+(1i*2*pi*freq(i)*C))';
    V_AC = (G+(1i*2*pi*freq(i)*C))\F;%F\(G+(1i*freq*C));%k*F 
    VoutAC(i)= abs(V_AC(5));
    
end

% AC case plot VO as a function of ω also plot the gain
figure(3)
plot(2*pi*freq,VoutAC)
Gain = db(VoutAC./F(6));
%xlim([0 150])
title('AC Output Voltage vs Frequency');
xlabel('Frequency (rad/s)')
ylabel('AC Output Voltage (V)')
figure(4)
plot(2*pi*freq,Gain)
%xlim([0 150])
title('Gain vs Frequency');
xlabel('Frequency (rad/s)')
ylabel('Gain (V/V)')


% This circuit is a 2nd-order circuit since there are two energy storing components. 
% Also, from the gain vs frequency response plot, it is a low pass filter.


Cv=zeros(100,1);
Cmatrix = zeros(7,7);
VoutCv = zeros(100,1);

w = pi;
 
for i = 1:100
     Cv(i) = 0.25 + 0.05*randn();
    Cmatrix(1,1) = Cv(i);
    Cmatrix(1,2) = -Cv(i);
    Cmatrix(2,6) = -L;
    Cmatrix(7,1) = -Cv(i);
    Cmatrix(7,2) = Cv(i);


    V_Cv = (G+(1i*pi*Cmatrix))\F;
    VoutCv(i)= abs(V_Cv(5));

end
% C case plot the gain as function of random perturbations on C using
% a normal distribution with std = .05 at ω = π. Do a histogram of the gain
GainCv = db(VoutCv./F(6));
figure(5)
histogram(GainCv,10);
title('Gain Distribution with different Capacitor values');
xlabel('Gain (V/V)')
ylabel('Number of Gains')


time = linspace(0,1,1000);
figure(6)
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

dt = 1e-3;
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

A = C/dt + G;
InvA = inv(A);

V1 = zeros(7,1);
V2 = zeros(7,1);
V3 = zeros(7,1);
F1 = zeros(7,1);
F2 = zeros(7,1);
F3 = zeros(7,1);


Vout1 = zeros(1000,1);
Vout2 = zeros(1000,1);
Vout3 = zeros(1000,1);
for i = 1:length(time)
    
    LastV1 = V1;
    LastV2 = V2;
    LastV3 = V3;
    
    F1(6) = Input1(i);
    F2(6) = Input2(i);
    F3(6) = Input3(i);

    V1 = A\(C*(LastV1/dt) + F1);
    Vout1(i,1) = V1(5,:);
    V2 = A\(C*(LastV2/dt) + F2);
    Vout2(i,1) = V2(5,:);
    V3 = A\(C*(LastV3/dt) + F3);
    Vout3(i,1) = V3(5,:);
end

figure(7)
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

frequency = linspace(-500,500,1000); % FFT Frequency Vector 

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

figure(8)
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

figure(9)
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

figure(10)
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

figure(11)
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