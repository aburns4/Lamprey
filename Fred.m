clear all
close all
clc
tic
n = 50;
m = 22;
alphaf = 4;
omega = 1;
omegaf = 1.08;
Aa = 8;
Ad = 6;
lambda_a = 20;
lambda_d = 20;
theta0 = [0:-1/n:-1+1/n 0];
dt=0.001;

options=odeset('RelTol',1e-4,'AbsTol',10e-7);

[T,Y] = ode45(@oscRHS,(0:dt:20),theta0,options,n,m,alphaf,omega,omegaf,Aa,Ad,lambda_a,lambda_d);
toc
plot(T,sin(2*pi*Y))
