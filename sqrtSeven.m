function x1 = sqrtSeven(x0)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
x0=3;
for i=0:3
    x1=1/2*(x0+7/x0);
    x0 = x1;
end

