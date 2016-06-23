pnr = 321;

A=dlmread('fort.38');

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
        area = area + segmentLength*(r(segNum,i) - l(segNum,i));
        areax = areax + segNum*(r(segNum,i) - l(segNum,i));
        areay = areay + 0.5*(r(segNum,i).^2 - l(segNum,i).^2);
    end
    xpoint(i) = areax/area;
    ypoint(i) = areay/area;
    fprintf('For time %d, the area is %0.3f\n',i,area);
    fprintf('For time %d, the point is (%0.2f,%0.2f)\n',i,xpoint(i),ypoint(i));
end


close all;
figure(1);

for i = 1:timeSteps
    clf('reset')
    
    plot(m(:,i))
    hold on
    plot(n(:,i))
    plot(l(:,i)) 
    plot(r(:,i))
    scatter(xpoint(i),ypoint(i));
    title(sprintf('Timestep is: %d',i));
    axis([0 350 5 20])
    pause(.1)
end