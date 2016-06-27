clear;
clc;
close all;

%BLUE CURVES
[xpos,ypos,xCenterMass,yCenterMass,xvel,yvel,timeSteps] = formatPositionFile('~/testForcingOtherSideStronger/fort.38');

%RED CURVES
%[xpos1,ypos1,xCenterMass1,yCenterMass1,xvel1,yvel1] = formatPositionFile('~/testForcing/fort.38');

for i = 1:timeSteps
    clf('reset')
    plot(xpos{1}(:,i),ypos{1}(:,i),'b')
    hold on
    plot(xpos{2}(:,i),ypos{2}(:,i),'b')
    plot(xpos{3}(:,i),ypos{3}(:,i),'b') 
    plot(xpos{4}(:,i),ypos{4}(:,i),'b')
    
%     plot(xpos1{1}(:,i),ypos1{1}(:,i),'r')
%     plot(xpos1{2}(:,i),ypos1{2}(:,i),'r')
%     plot(xpos1{3}(:,i),ypos1{3}(:,i),'r')
%     plot(xpos1{4}(:,i),ypos1{4}(:,i),'r')
    
    axis([-20 20 0 12])
    scatter(xCenterMass(i),yCenterMass(i),'b');
    quiver(xCenterMass(i),yCenterMass(i),xvel(i),yvel(i),'Color','b','MaxHeadSize',0.15,'AlignVertexCenters','on','LineWidth',1.0);
    
%     scatter(xCenterMass1(i),yCenterMass1(i),'r');
%     quiver(xCenterMass1(i),yCenterMass1(i),xvel1(i),yvel1(i),'Color','r','MaxHeadSize',0.15,'AlignVertexCenters','on','LineWidth',1.0);
    
    title(sprintf('Timestep is: %d',i));
    pause(.01)
end
plot(xCenterMass,yCenterMass,'b');
% plot(xCenterMass1,yCenterMass1,'r');