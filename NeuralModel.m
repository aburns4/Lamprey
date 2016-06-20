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

Aa = 1;
Ad = 1;
lambda_a = (-1*log(0.05))^(-1);
lambda_d = (-1*log(0.02))^(-1);
omega_f = 1; %Forcing frequency

n = 1; %Number of oscillators
m = 0; %Forcing position
dt = 0.001; %Time step

alpha_f = 5; %Forcing strength

v_0 = zeros(6*n+1,1); %Initial Conditions
%v_0 = (2*rand(6*n+1,1))-1;

options=odeset('RelTol',1e-4,'AbsTol',1e-7);

%alpha_r = CouplingFunction(n,m,Aa,Ad,lambda_a,lambda_d,sigma);
%alpha_r is the strength of the connections for different r
%rows are of r and columns are [L to C, E to C, other]

alpha_r = zeros(2*n-1,3);
for r = 1:(2*n-1)
    k = r-n;
    alpha_r(r,:)=str(Aa,Ad,lambda_a,lambda_d,k);
end
    

[T,Y] = ode45(@neuralFunc,(0:dt:10),v_0,options,n,m,G_R,G_T,G_0,V_syn,G_f,V_synec,sigma,alpha_f,alpha_r,omega_f);
[T,Y] = ode45(@neuralFunc,(0:dt:3),Y(end,:),options,n,m,G_R,G_T,G_0,V_syn,G_f,V_synec,sigma,alpha_f,alpha_r,omega_f);

% plot(T,Y(:,1),'b');
% hold on
% plot(T,Y(:,2),'color',[0 0.5 0]);
% plot(T,Y(:,3),'r');
% plot(T,Y(:,4),'b--');
% plot(T,Y(:,5),'--','color',[0 0.5 0]);
% plot(T,Y(:,6),'r--');
% hold off
% title('Neural Model Simulation','FontName', 'Times New Roman',  'FontSize', 20);
% legend('Left E','Left L','Left C', 'Right E','Right L','Right C', 'FontName', 'Times New Roman', 'FontSize', 20);
% xlabel('Time in Seconds', 'FontName', 'Times New Roman', 'FontSize', 20);
% ylabel('Cell Voltage V_{ij}', 'FontName', 'Times New Roman', 'FontSize', 20);
%axis([0 2.7 -1 1]);
figure(1)
plot(T,Y(:,1:6*n));

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

%CALCULATES THE PERIOD AND NATURAL FREQUENCY OF THE SYSTEM
E_zeros=find(Y(1:end-1,1)<0 & Y(2:end,1)>=0);
L_zeros=find(Y(1:end-1,2)<0 & Y(2:end,2)>=0);
C_zeros=find(Y(1:end-1,3)<0 & Y(2:end,3)>=0);
y1=E_zeros(2)-E_zeros(1);
y2=L_zeros(2)-L_zeros(1);
y3=C_zeros(2)-C_zeros(1);
period = (y1+y2+y3)*dt/3;
freq = 1/period;
fprintf('The period is %f and the frequency is %f \n',period,freq);


toc