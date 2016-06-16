function a = evenOdd(x)
%if the input is even it halves it, otherwise it is odd and it is doubled

if(mod(x,2) == 0);
    a = x/2;
else
    a = x*2;

end

