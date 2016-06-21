%run forwardEulerAgain

close all;
clear;
clc;

t_0=0;
t_f=76;
dt = 0.0005; %Time step

G_R = 3.5; %Resting Conductance
G_T = [0.875, 0.35, 3.5]; %Tonic Excitatory Conductance [E,L,C]
G_0 = [15, 35]; %Maximal Synaptic Conductance of Intersegmental Connection
                %[L to C, other]
V_syn = [1, -1]; %Synaptic reversal potential for intersegmental connection
                 %[excitatory, inhibatory]
G_f = 1; %Maximal synaptic conductance of EC connections
V_synec = [1, -1]; %Synaptic reversal potential for EC connections [excitatory, inhibatory]
sigma = 0.05; %Parameter for smooth threshold function

omega_f = 1; %Forcing frequency

n = 5; %Number of oscillators
m = 0; %Forcing position

alpha_f = 5; %Forcing strength

v_0 = zeros(6*n+1,1); %Initial Conditions
%v_0 = (2*rand(6*n+1,1))-1;

alpha_r=[.05 1 .02];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[T,Y] = forwardEulerAgain(t_0,t_f,dt,v_0,n,m,G_R,G_T,G_0,V_syn,G_f,V_synec,sigma,alpha_f,alpha_r,omega_f);

%plot(T,Y(:,1:6*n))

figure(1)
for i= 1:size(Y,2)
    if mod(i,6)==0
        plot(T,Y(:,i))
        hold on
    end
end

figure(2)
for i= 1:size(Y,2)-6
    if mod(i,6)==0
        plot(T,Y(:,i+1))
        hold on
    end
end

figure(3)
for i= 1:size(Y,2)-6
    if mod(i,6)==0
        plot(T,Y(:,i+2))
        hold on
    end
end


%CALCULATES THE PERIOD AND NATURAL FREQUENCY OF THE SYSTEM
period = zeros(1,n);
j=0;
for i = 1:6:6*n
    j=j+1;
        
        E_zeros=find(Y(1:end-1,i)<0 & Y(2:end,i)>=0);
        L_zeros=find(Y(1:end-1,i+1)<0 & Y(2:end,i+1)>=0);
        C_zeros=find(Y(1:end-1,i+2)<0 & Y(2:end,i+2)>=0);
        y1=E_zeros(end)-E_zeros(end-1);
        y2=L_zeros(end)-L_zeros(end-1);
        y3=C_zeros(end)-C_zeros(end-1);
        period(j) = (y1+y2+y3)*dt/3;
        freq(j) = 1/period(j);
end






%Plots Phase Lag
% figure(2)
% for i = 1:(2*n)-8
%     E1 = Y(:,i); %first E cell, left side
%     L1 = Y(:,i+1);%first L cell, left side
%     C1 = Y(:,i+2);%first C cell, left side
%     E2 = Y(:,i+6);%second E cell, left side
%     L2 = Y(:,i+7);%second L cell, left side
%     C2 = Y(:,i+8);%second C cell, left side
%     
%     [A1, IA1] = max(E1);
%     [A2, IA2] = max(E2);
%     
%     [B1, IB1] = max(L1);
%     [B2, IB2] = max(L2);
%     
%     [D1, ID1] = max(C1);
%     [D2, ID2] = max(C2);
%     
%     Ediff = (IA2 - IA1)*dt;
%     Ldiff = (IB2 - IB1)*dt;
%     Cdiff = (ID2 - ID1)*dt;
%     scatter(T(IA2),Ediff);
%     hold on
%     scatter(T(IB2),Ldiff);
%     
%     scatter(T(ID2),Cdiff);
%     
% end