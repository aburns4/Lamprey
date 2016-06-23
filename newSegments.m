function [final_E] = newSegments(eCell,period,dt,ratio)

active = find(eCell == 1); %finds the indices of all 1's, I want to know how long cell is active  

duty_cycle = (size(active,1)*dt)/period; %calculate duty cycle
base_lag = duty_cycle/(ratio-1); %calculate a base phase lag dependent on duty cycle and how many new segments there are
new_E = zeros(size(eCell,1)+floor((base_lag*(ratio-1))/dt),1); %where I will store new E cell values

%matrix of all new E cells per old segment (right now its 4 new cells per old
%E cell
final_E = zeros(size(eCell,1)+floor((base_lag*(ratio-1))/dt),ratio); 


%looping through every new segment to prescribe phase lag
for k = 1:ratio
    %how many time steps for the prescribed phase lag (make everything
    %before this value zero and then set elements to old eCell
    
    time_steps = floor((base_lag*(k-1))/dt); %first column no phase lag, increases with k
    
    %%%%%%THE PROBLEM IS MOST LIKELY IN THE LINE BELOW%%%%%%%%%%%
    %%%IT IS NOT SHIFTING WHERE IT IS SETTING IT TO ECELL EVEN THOUGH
    %%%TIMESTEP DOES CHANGE SO WHAT THE FUCK?!%%%%%%%%%%%%%%%%%
    new_E(time_steps+1:time_steps+size(eCell,1))= eCell;
    
    %saves all new segments (per that one old e cell) in a 2d matrix 
    final_E(:,k) = new_E;
end


end

