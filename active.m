function thetaValues = active(t0,tf,theta0,dt,n,dc,w,fposition,wf,thetaf0,alphaf)
%this calls the Euler method on our input values of active
    thetaValues = Euler(t0,tf,theta0,dt,n,w,fposition,wf,thetaf0,alphaf);

    %this loops through every value of the matrix and checks if it is above
    %the necessary value to have the proper equal dc per point
    for i = 1:n
        for j = 1:size(thetaValues,2)
            if thetaValues(i,j) >= sin(2*pi*w*(.25*w+dc/2)) 
                %if true, replace with one as "active"
                thetaValues(i,j) = 1;
            else
                %replace with zero as "inactive"
                thetaValues(i,j) = 0;
            end
        end
    end
end

