clear;
clc;
close all;

pnr = 321;

A=dlmread('~/testForcing1/fort.38');

timeSteps=size(A,1)/(4*pnr);

X=reshape(A(:,1),size(A,1)/timeSteps,timeSteps);
Y=reshape(A(:,2),size(A,1)/timeSteps,timeSteps);

Xm=X(1:pnr,:);          Ym=Y(1:pnr,:);
Xn=X(1+pnr:2*pnr,:);    Yn=Y(1+pnr:2*pnr,:);
Xl=X(1+2*pnr:3*pnr,:);  Yl=Y(1+2*pnr:3*pnr,:);
Xr=X(1+3*pnr:end,:);    Yr=Y(1+3*pnr:end,:);

m=zeros(pnr,timeSteps);
n=zeros(pnr,timeSteps);
l=zeros(pnr,timeSteps);
r=zeros(pnr,timeSteps);

for i = 1:timeSteps
    m(:,i)=sqrt(Xm(:,i).^2 + Ym(:,i).^2);
    n(:,i)=sqrt(Xn(:,i).^2 + Yn(:,i).^2);
    l(:,i)=sqrt(Xl(:,i).^2 + Yl(:,i).^2);
    r(:,i)=sqrt(Xr(:,i).^2 + Yr(:,i).^2);
end

segmentLength = 1;

xpoint = zeros(1,timeSteps);
ypoint = zeros(1,timeSteps);

for i = 1:timeSteps
    area = 0;
    areax = 0;
    areay = 0;
    for segNum = 1:pnr
        area = area + segmentLength*(Yr(segNum,i) - Yl(segNum,i));
        areax = areax + Xm(segNum,i)*(Yr(segNum,i) - Yl(segNum,i));
        areay = areay + 0.5*(Yr(segNum,i).^2 - Yl(segNum,i).^2);
    end
    xpoint(i) = areax/area;
    ypoint(i) = areay/area;
    fprintf('For time %d, the area is %0.3f\n',i,area);
    fprintf('For time %d, the point is (%0.2f,%0.2f)\n',i,xpoint(i),ypoint(i));
end

figure(1);
dt = 0.025;
xvel = zeros(1,timeSteps);
yvel = zeros(1,timeSteps);

for i = 1:timeSteps
    clf('reset')
    
    plot(Xm(:,i),Ym(:,i))
    hold on
    plot(Xn(:,i),Yn(:,i))
    plot(Xl(:,i),Yl(:,i)) 
    plot(Xr(:,i),Yr(:,i))
    axis([-20 20 0 12])
    scatter(xpoint(i),ypoint(i));
    if(i>1 && i<timeSteps)
        xvel(i) = (xpoint(i+1)-xpoint(i-1))/(2*dt);
        yvel(i) = (ypoint(i+1)-ypoint(i-1))/(2*dt);
        quiver(xpoint(i),ypoint(i),xvel(i),yvel(i),'MaxHeadSize',0.15,'AlignVertexCenters','on','LineWidth',1.0);
    elseif(i == timeSteps)
        quiver(xpoint(i-1),ypoint(i-1),xvel(i-1),yvel(i-1),'MaxHeadSize',0.15,'AlignVertexCenter','on','LineWidth',1.0);
    end
    title(sprintf('Timestep is: %d',i));
    pause(.1)
end
plot(xpoint,ypoint);