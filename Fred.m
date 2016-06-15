
%clear all
%close all
clc
tic
n = 10;
m = 0;
alphaf = 4;
omega = 0.5;
omegaf = 1;
Aa = 1;
Ad = 10;
lambda_a = 40;
lambda_d = 5;
%theta0 = [0:-1/n:-1+1/n -(m-1)/n];
theta0 = zeros(n+1,1);
dt=0.001;

options=odeset('RelTol',1e-4,'AbsTol',10e-7);

[H_r, H_forcing, Hs1]=CouplingFunction(n,m,Aa,Ad,lambda_a,lambda_d);

[T,Y] = ode45(@oscRHS,(0:dt:5),theta0,options,n,m,alphaf,omega,omegaf,Aa,Ad,lambda_a,lambda_d,H_r,H_forcing,Hs1);
toc

close all;

plot(T,sin(2*pi*Y(:,1:n)));
%plot(T,Y(:,1:n));

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
