function [fs,theta,STRENGTH] = Euler(t0,tf,theta0,dt,n,w,fposition,wf,thetaf0,alphaf)
tic
%Sets initial time to t0
t(1)=t0;

%Creates initial conditions with phase lag and places it into theta
if(size(theta0,1)==1 && size(theta0,2)==1)
    initial=zeros(n,1);
    for i=1:n
        initial(i) = theta0 - (i-1)/n;
    end
    theta(:,1)=initial;
elseif(size(theta0,1)==1 && size(theta0,2)==n)
    theta(:,1)=theta0';
elseif(size(theta0,1)==n && size(theta0,2)==1)
    theta(:,1)=theta0;
end

%Initialize the forcing term
counter=1;
forcingTheta=zeros(n,1);
forcingOmega=zeros(n,1);
for i=1:n
    if(i>fposition(size(fposition,2)))
        break;
    elseif(i==fposition(counter))
        forcingTheta(i)=thetaf0(counter);
        forcingOmega(i)=wf(counter);
        counter=counter+1;
    end
end

%Finds number of steps necessary
N=floor((tf-t0)/dt);
t=zeros(1,N-1);
for i=1:(N-1)
    %Use allToAll function for all to all coupling
    %Use nearestNeighbor function for nearest neighbor coupling
    [temp,STRENGTH]=allToall(t(i),theta(:,i),n,w,forcingTheta(:,i),forcingOmega,alphaf);
    forcingTheta(:,i+1)=forcingTheta(:,i)+dt*forcingOmega;
%    %Plot the surface with connection distance
%     if(mod(i,30)==0)
%         surf(STRENGTH);
%         xlabel('Cell');
%         ylabel('Cell connection');
%         zlabel('Cell distance');
%         title(sprintf('Time = %.3f',i*dt));
%         pause(0.001);
%     end
    
    theta(:,i+1) = theta(:,i) + dt*(temp);
    t(i+1) = t(i) + dt;
end
%Transforms the theta matrix into a signal matrix
fs=sin(2*pi*theta);
toc
end
