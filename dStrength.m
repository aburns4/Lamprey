function alpha = dStrength(Ad,ld,k)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    alpha = Ad*exp(-abs(k)/ld);
end

