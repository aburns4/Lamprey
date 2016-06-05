%this is a function to compute forcing strength, aplha_{ij}
function aij = fStrength(i, j, n, thetaj, thetai)
    %here we calculate phase lag as the difference in chain cell over total
    %number of cells
    lag = (i-j)/n;
    
    %if i-j is less than zero it is ascending
    if (i - j) < 0
        alpha = str(12,1/log(12/10.1),abs(i-j));
        
    %same cell if i==j, so do nothing
    elseif i == j
        alpha = 0;
        
    %else we have i-j > 0 which is descending
    else
        alpha = str(12,1/log(12/10.0),abs(i-j));
    
    end
        
    %calculate a_ij with given sine equation
    aij = alpha*sin(2*pi*(thetaj-thetai-lag));
 
end

