tic

close all;
clear;
clc;

G_R = 3.5; %Resting Conductance
G_T = [0.875, 0.35, 3.5]; %Tonic Excitatory Conductance [E,L,C]
G_0 = [15, 35]; %Maximal Synaptic Conductance of Intersegmental Connection
                %[L to C, other]
V_syn = [1, -1]; %Synaptic reversal potential for intersegmental connection
                 %[excitatory, inhibatory]
G_f = 1; %Maximal synaptic conductance of EC connections
V_synec = [1, -1]; %Synaptic reversal potential for EC connections [excitatory, inhibatory]
sigma = 0.05; %Parameter for smooth threshold function

n = 2;
m = 0;
alphaf = 1;
omega = 1;
omegaf = 1.0005;
Aa = 1;
Ad = 1;
lambda_a = (-1*log(0.05))^(-1);
lambda_d = (-1*log(0.02))^(-1);
v_0 = zeros(6*n+1,1); %Initial Conditions
dt=0.001;

options=odeset('RelTol',1e-4,'AbsTol',10e-7);

[T,Y] = ode45(@neuralFunc,(0:dt:50),v_0,options,n,m,G_R,G_T,G_0,V_syn,G_f,V_synec,sigma,alpha_f,alpha_r,omega_f);

v_0 = Y(end,:);
[T,Y] = ode45(@neuralFunc,(0:dt:3),Y(end,:),options,n,m,G_R,G_T,G_0,V_syn,G_f,V_synec,sigma,alpha_f,alpha_r,omega_f);

%CALCULATES THE PERIOD AND NATURAL FREQUENCY OF THE SYSTEM
E_zeros=find(Y(1:end-1,1)<0 & Y(2:end,1)>=0);
L_zeros=find(Y(1:end-1,2)<0 & Y(2:end,2)>=0);
C_zeros=find(Y(1:end-1,3)<0 & Y(2:end,3)>=0);
y1=E_zeros(2)-E_zeros(1);
y2=L_zeros(2)-L_zeros(1);
y3=C_zeros(2)-C_zeros(1);
periodNatural = (y1+y2+y3)*dt/3;
freqNatural = 1/periodNatural;
fprintf('The period is %f and the frequency is %f \n',periodNatural,freq);


for forcingPosition=1:n
    for omega_f= freqNatural:.0001:2
        [T,Y] = ode45(@neuralFunc,(0:dt:50),v_0,options,n,forcingPosition,G_R,G_T,G_0,V_syn,G_f,V_synec,sigma,alpha_f,alpha_r,omega_f);
        [T,Y] = ode45(@neuralFunc,(0:dt:3),Y(end,:),options,n,forcingPosition,G_R,G_T,G_0,V_syn,G_f,V_synec,sigma,alpha_f,alpha_r,omega_f);

        E_zeros=find(Y(1:end-1,1)<0 & Y(2:end,1)>=0);
        L_zeros=find(Y(1:end-1,2)<0 & Y(2:end,2)>=0);
        C_zeros=find(Y(1:end-1,3)<0 & Y(2:end,3)>=0);
        y1=E_zeros(2)-E_zeros(1);
        y2=L_zeros(2)-L_zeros(1);
        y3=C_zeros(2)-C_zeros(1);
        periodForcing = (y1+y2+y3)*dt/3;
        freqAfterForcing = 1/period;
        
        if (abs(freqAfterForcing-omega_f) >= .000001)
            etrRanges(forcingPosition) = abs(omega_f - freqNatural);
            break
        end
    end
end

toc