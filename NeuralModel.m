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
omega_f=1; %Forcing frequency

n = 10; %Number of oscillators
m = 0; %Forcing position
dt = 0.001; %Time step

alpha_f = 0; %Forcing strength

%v_0 = zeros(6*n+1,1); %Initial Conditions
v_0 = (2*rand(6*n+1,1))-1;

options=odeset('RelTol',1e-4,'AbsTol',10e-7);

alpha_r = CouplingFunction(n,m,Aa,Ad,lambda_a,lambda_d);
%alpha_r is the strength of the connections for different r
%rows are of r and columns are [L to C, E to C, other]

[T,Y] = ode45(@neuralFunc,(0:dt:2),v_0,options,n,m,G_R,G_T,G_0,V_syn,G_f,V_synec,sigma,alpha_f,alpha_r,omega_f);





