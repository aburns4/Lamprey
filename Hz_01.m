function [phi0, dH] = Hz_01(phi, H)
% [phi0 dH] = Hz_01(phi, H)
%
% Finds zeros of coupling functions.

m    = size(H, 2);
phi0 = zeros(m, 1);
dH   = zeros(m, 1);


for i = 1 : m,
 j = find(H(1:end-1,i) < 0 & H(2:end,i) >= 0); %look for zero crossings since only have grid values for H

 if isempty(j)
   phi0(i) = NaN;
 else
   phi0(i) = phi(j) - (phi(j+1) - phi(j))*H(j,i)/(H(j+1,i) - H(j,i)); %like Newton step
   dH(i)   = (H(j+1,i) - H(j,i))/(phi(j+1) - phi(j)); %derviative
 end
end
