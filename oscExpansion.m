function [final_activation] = oscExpansion(new_n)

load activationLevels.mat; %load activation matrix and n

old_n = n; %how many segments there were- probably around 20, where the limit was reached

ratio = floor(new_n/old_n); %how many new E cells correlate to an old E cell

new_EAll = cell(1,old_n); %holds the subcells of matrices, each being a new set of E values

final_activation = zeros(size(activationLR,1),2*new_n); %where all new activation levels will be store

for i = 1:size(activationLR,2)-1
    %%%%%%%%%%%%%Call function to compute more activation levels%%%%%%%%%%%%%
    new_E = newSegments(activationLR(:,i),period,dt,ratio);%a matrix of E cells corresponding to one old E cell, # of columns = ratio above
    new_EAll{i} = new_E; %store matrix as subcell 
end

for i = 1:old_n-1
    for k = 1:ratio
        final_activation(:,i) = new_EAll{i}(:,k); %left E cell
        final_activation(:,i+1) = new_EAll{i+1}(:,k); %right E cell
    end
end

end