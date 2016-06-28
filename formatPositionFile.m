function [xpos,ypos,xCenterMass,yCenterMass,xvel,yvel,timeSteps] = formatPositionFile(fileName)
    %xpos and ypos are made of 4 cells:
    %       cell 1 -> Middle (m)
    %       cell 2 -> Middle (n)
    %       cell 3 -> Left   (l)
    %       cell 4 -> Right  (r)

    A=dlmread(fileName);
    pnr=321;
    timeSteps=size(A,1)/(4*pnr);

    X=reshape(A(:,1),size(A,1)/timeSteps,timeSteps);
    Y=reshape(A(:,2),size(A,1)/timeSteps,timeSteps);

    xpos{1}=X(1:pnr,:);          ypos{1}=Y(1:pnr,:);
    xpos{2}=X(1+pnr:2*pnr,:);    ypos{2}=Y(1+pnr:2*pnr,:);
    xpos{3}=X(1+2*pnr:3*pnr,:);  ypos{3}=Y(1+2*pnr:3*pnr,:);
    xpos{4}=X(1+3*pnr:end,:);    ypos{4}=Y(1+3*pnr:end,:);
    
    segmentLength = 1;

    xCenterMass = zeros(1,timeSteps);
    yCenterMass = zeros(1,timeSteps);
    
    for i = 1:timeSteps
        area = 0;
        areax = 0;
        areay = 0;
        for segNum = 1:pnr
            area = area + segmentLength*(ypos{4}(segNum,i) - ypos{3}(segNum,i));
            areax = areax + xpos{1}(segNum,i)*(ypos{4}(segNum,i) - ypos{3}(segNum,i));
            areay = areay + 0.5*(ypos{4}(segNum,i).^2 - ypos{3}(segNum,i).^2);
        end
        xCenterMass(i) = areax/area;
        yCenterMass(i) = areay/area;
        %fprintf('For time %d, the area is %0.3f\n',i,area);
        %fprintf('For time %d, the point is (%0.2f,%0.2f)\n',i,xCenterMass(i),yCenterMass(i));
    end
    
    xvel = zeros(1,timeSteps);
    yvel = zeros(1,timeSteps);
    dt = 0.025;
    
    for i = 1:timeSteps
        if(i>1 && i<timeSteps)
            xvel(i) = (xCenterMass(i+1)-xCenterMass(i-1))/(2*dt);
            yvel(i) = (yCenterMass(i+1)-yCenterMass(i-1))/(2*dt);
        end
    end
    xvel(1) = xvel(2);
    xvel(end) = xvel(end-1);
end