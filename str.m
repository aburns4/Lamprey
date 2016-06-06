function alpha = str(A,l,k)
%generalized function for ascending OR descending coupling strength
    alpha = A*exp(-abs(k)/l);
end

