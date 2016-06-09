tic
n = 10;
m = 1;
alphaf = 0.01;
omega = 1;
omegaf = 1.0005;
Aa = 0.0004;
Ad = 0.0002;
lambda_a = 4;
lambda_d = 4;
theta0 = [0:-1/n:-1+1/n 0];
dt=0.001;

options=odeset('RelTol',1e-4,'AbsTol',10e-7);

%[T,Y] = ode45(@oscRHS,(0:dt:35),theta0,options,n,m,alphaf,omega,omegaf,Aa,Ad,lambda_a,lambda_d);

segment=1:n;
for forcingPos=1:n
    for omegaf=1:0.0001:1.0006
        [T,Y]=ode45(@oscRHS,(0:dt:30),theta0,options,n,forcingPos,alphaf,omega,omegaf,Aa,Ad,lambda_a,lambda_d,H,shift);
        for k=1:20
            A=(Y(size(Y,1)-k+1,:)-Y(size(Y,1)-k-1,:))/(2*dt);
            B(k)=mean(A);
        end
        avgFreq=mean(B);
        clear T Y B;
        avgFreq=round(avgFreq,4);
        if(avgFreq-omegaf~=0)
            freq(forcingPos)=omegaf-omega;
            break;
        end
    end
end
toc
%plot(segment,freq,'b',segment,-1*freq,'b');
        