function alpha = aStrength(Aa,la,k)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    alpha = Aa*exp(-abs(k)/la);
end

