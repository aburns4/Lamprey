%run forwardEulerAgain

close all;
clear;
clc;

t_0=0;
t_f=10;
dt = 0.001; %Time step

G_R = 3.5; %Resting Conductance
G_T = [0.875, 0.35*4.05, 3.5]; %Tonic Excitatory Conductance [E,L,C]
G_0 = [15, 35]; %Maximal Synaptic Conductance of Intersegmental Connection
                %[L to C, other]
V_syn = [1, -1]; %Synaptic reversal potential for intersegmental connection
                 %[excitatory, inhibatory]
G_f = 1; %Maximal synaptic conductance of EC connections
V_synec = [1, -1]; %Synaptic reversal potential for EC connections [excitatory, inhibatory]
sigma = 0.05; %Parameter for smooth threshold function

omega_f = 1; %Forcing frequency

n = 15; %Number of oscillators
m = 0; %Forcing position

alpha_f = 5; %Forcing strength

v_0 = zeros(6*n+1,1); %Initial Conditions
%v_0 = (2*rand(6*n+1,1))-1;

alpha_r=[.05 1 .02];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[T,Y] = forwardEulerAgain(0,30,dt,v_0,n,m,G_R,G_T,G_0,V_syn,G_f,V_synec,sigma,alpha_f,alpha_r,omega_f);
[T,Y] = forwardEulerAgain(t_0,t_f,dt,Y(end,:),n,m,G_R,G_T,G_0,V_syn,G_f,V_synec,sigma,alpha_f,alpha_r,omega_f);


LeftECellMat = zeros(size(Y,1),n);
RightECellMat = zeros(size(Y,1),n);

for i=1:n
    LeftECellMat(:,i) = Y(:,6*(i-1)+1);
    RightECellMat(:,i) = Y(:,6*(i-1)+4);
    for j=1:size(Y,1)
        if(LeftECellMat(j,i)>0)
            LeftECellMat(j,i)=1;
        else
            LeftECellMat(j,i)=.2;
        end
        if(RightECellMat(j,i)>0)
            RightECellMat(j,i)=1;
        else
            RightECellMat(j,i)=.2;
        end
    end
end

% figure(1);
% plot(LeftECellMat);
% figure(2);
% plot(RightECellMat);


for t=1:length(T)
    if mod(t,30)==0
        figure(1)
       
        plot(LeftECellMat(t,:))
        hold on
        plot(-RightECellMat(t,:))
        
        axis([.9 n+.1 -1.1 1.1])
        title(sprintf('Time = %.2f',t*dt));
        pause(.09)
        
        clf('reset')
    end
end