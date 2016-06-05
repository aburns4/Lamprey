%this defines our function for the Euler method for neighbor coupling

function [f,FORCING] = neighbor(t,theta,n,w)

for i=1:n   
    %case: first cell in chain, only 1 connection
     if i ==1
        f(i) = w + fStrength(i,i+1,n,theta(i+1),theta(i));
        FORCING(i,i+1) = fStrength(i,i+1,n,theta(i+1),theta(i));
    %case: end cell in chain, only 1 connection
    elseif i == n
        f(i) = w + fStrength(i,i-1,n,theta(i-1),theta(i));
        FORCING(i,i-1) = fStrength(i,i-1,n,theta(i-1),theta(i));
        
     %case: middle cell in chain, so it has 2 connections
     else   
         f(i) = w + fStrength(i,i+1,n,theta(i+1),theta(i)) + fStrength(i,i-1,n,theta(i-1),theta(i));  
         FORCING(i,i+1) = fStrength(i,i+1,n,theta(i+1),theta(i));
         FORCING(i,i-1) = fStrength(i,i-1,n,theta(i-1),theta(i));
     end
     
    end
 f = f';
end

