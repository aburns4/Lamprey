function [] = plotActive(t0,tf,theta0,dt,n,dc,w,fposition,wf,thetaf0,alphaf)

    thetaValues = active(t0,tf,theta0,dt,n,dc,w,fposition,wf,thetaf0,alphaf);
    %this goes through the columns which are different time stamps and then
    %plots a scatterplot of each row's values for active
    for i=1:size(thetaValues,2)
        if (mod(i,30)==0)
        y = thetaValues(:,i);
        scatter(1:n,y);
        line(1:n,y);
        xlabel('section');
        ylabel('activation');
        title(sprintf('Time = %.2f',i*dt));
        pause(0.01);
        end
    end 
end

