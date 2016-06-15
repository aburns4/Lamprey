
%clear all
%close all
clc
tic
n = 5;
m = 3;
alphaf = 0.02;
omega = 1;
omegaf = 1;
Aa = 0.006;
Ad = 0.0004;
lambda_a = 0.75;
lambda_d = 4;
theta0 = [0:-1/n:-1+1/n -(m-1)/n];
dt=0.001;

options=odeset('RelTol',1e-4,'AbsTol',10e-7);

[H_r, H_forcing]=CouplingFunction(n,m,Aa,Ad,lambda_a,lambda_d);

[T,Y] = ode45(@oscRHS,(0:dt:10),theta0,options,n,m,alphaf,omega,omegaf,Aa,Ad,lambda_a,lambda_d,H_r,H_forcing);
toc

close all;

%plot(T,sin(2*pi*Y(:,1:n)));
plot(T,Y(:,1:n));

%plots  voltage with forcing. To change plotting speed change mod values 
% for t = 1:(20/dt)
%     
%     if mod(t,55)==0
%      plot(1:n,sin(2*pi*Y(t,1:n)))
%         z=sin(2*pi*omegaf*t*dt);
%     
%          hold on
%     
%      scatter(m,z)
%      
%       xlabel('section');
%         ylabel('Voltage');
%         title(sprintf('Time = %.2f',t*dt));
%      
%         pause(dt)
%  
%         clf('reset')
%     end
% end
