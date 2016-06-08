tic
n = 50;
m = 1;
alphaf = 4;
omega = 1;
omegaf = 1;
Aa = 8;
Ad = 6;
lambda_a = 20;
lambda_d = 20;
theta0 = [0:-1/n:-1+1/n 0];
dt=0.001;

options=odeset('RelTol',1e-4,'AbsTol',10e-7);

%[T,Y] = ode45(@oscRHS,(0:dt:35),theta0,options,n,m,alphaf,omega,omegaf,Aa,Ad,lambda_a,lambda_d);

segment=1:n;
for forcingPos=1:n
    for omegaf=1:0.01:1.15
        [T,Y]=ode45(@oscRHS,(0:dt:60),theta0,options,n,forcingPos,alphaf,omega,omegaf,Aa,Ad,lambda_a,lambda_d);
        for k=1:20
            A=(Y(size(Y,1)-k+1,:)-Y(size(Y,1)-k-1,:))/(2*dt);
            B(k)=mean(A);
        end
        avgFreq=mean(B);
        clear T Y B;
        avgFreq=round(avgFreq,2);
        if(avgFreq-omegaf~=0)
            freq(forcingPos)=omegaf-omega;
            break;
        end
    end
end
toc
%plot(segment,freq,'b',segment,-1*freq,'b');
        