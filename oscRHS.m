function dy = oscRHS(t,theta,n,m,alphaf,omega,omegaf,Aa,Ad,lambda_a,lambda_d,H,shift)
    dy = zeros(n+1,1);
    
    for i = 1:n
        couplingstr = 0;
        for j = 1:n
            k = i - j;
            if k < 0
                alpha_k = Aa*exp(-abs(k)/lambda_a);
            elseif k > 0
                alpha_k = Ad*exp(-abs(k)/lambda_d);
            else
                alpha_k = 0;
            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             phaseDiff=theta(j)-theta(i);%-(i-j)/n;
             if(phaseDiff<0)
                 phaseDiff=mod(abs(phaseDiff),pi)*-1;
             else
                 phaseDiff=mod(phaseDiff,pi);
             end
            couplingfctn=interp1(-shift{4},H{4},phaseDiff,'spline');
            couplingstr = couplingstr + alpha_k.*(couplingfctn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %couplingstr=couplingstr+alpha_k.*sin(theta(j)-theta(i)-(i-j)/n);
        end
        if i == m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            phaseDiff=theta(n+1)-theta(i);
            if(phaseDiff<0)
                phaseDiff=mod(abs(phaseDiff),pi)*-1;
            else
                phaseDiff=mod(phaseDiff,pi);
            end
            couplingfctn=interp1(-shift{4},H{4},phaseDiff,'spline');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %couplingstr = couplingstr + alphaf*sin(2*pi*(theta(n+1) - theta(i)));
            couplingstr=couplingstr+alphaf*(couplingfctn);
        end
        dy(i) = omega + couplingstr; 
    end
    dy(n+1) = omegaf;
end

