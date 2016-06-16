G_R = 3.5; %Resting Conductance
G_T = [0.875, 0.35, 3.5]; %Tonic Excitatory Conductance [E,L,C]
G_0 = [15, 35]; %Maximal Synaptic Conductance of Intersegmental Connection
                %[L to C, other]
V_syn = [1, -1]; %Synaptic reversal potential for intersegmental connection
                 %[excitatory, inhibatory]
G_f = 1; %Maximal synaptic conductance of EC connections
V_synec = [1, -1]; %Synaptic reversal potential for EC connections [excitatory, inhibatory]
sigma = 0.05; %Parameter for smooth threshold function

Aa = 0.0004;
Ad = 0.0002;
lambda_a = 4;
lambda_d = 4;
omega_f=.7407; %Forcing frequency

n = 1; %Number of oscillators
m = 1; %Forcing position
dt = 0.001; %Time step

alpha_f = 0.02; %Forcing strength

%v_0 = zeros(6*n+1,1); %Initial Conditions
%v_0 = (2*rand(6*n+1,1))-1;
v_0 = [0.2 0.5 0.19 -0.55 -0.6 -0.3 0];

options=odeset('RelTol',1e-4,'AbsTol',1e-7);

alpha_r = CouplingFunction(n,m,Aa,Ad,lambda_a,lambda_d);
%alpha_r is the strength of the connections for different r
%rows are of r and columns are [L to C, E to C, other]

[T,Y] = ode45(@neuralFunc,(0:dt:5),v_0,options,n,m,G_R,G_T,G_0,V_syn,G_f,V_synec,sigma,alpha_f,alpha_r,omega_f);

plot(T,Y(:,1),'b');
hold on
plot(T,Y(:,2),'color',[0 0.5 0]);
plot(T,Y(:,3),'r');
plot(T,Y(:,4),'b--');
plot(T,Y(:,5),'--','color',[0 0.5 0]);
plot(T,Y(:,6),'r--');
hold off
title('Neural Model Simulation','FontName', 'Times New Roman',  'FontSize', 20);
legend('Left E','Left L','Left C', 'Right E','Right L','Right C', 'FontName', 'Times New Roman', 'FontSize', 20);
xlabel('Time in Seconds', 'FontName', 'Times New Roman', 'FontSize', 20);
ylabel('Cell Voltage V_{ij}', 'FontName', 'Times New Roman', 'FontSize', 20);
axis([0 2.7 -1 1]);


