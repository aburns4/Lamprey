function dy = oscRHS(t,theta,n,m,alphaf,omega,omegaf,Aa,Ad,lambda_a,lambda_d,H_r,H_forcing,Hs1)

%load CouplingFunction.mat;

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
        phaseDiff=theta(j)-theta(i);
        if(phaseDiff<0)
            phaseDiff=mod(abs(phaseDiff),0.5)*-1;
        else
            phaseDiff=mod(phaseDiff,0.5);
        end
        couplingfctn = interp1(-0.5:0.01:0.5,Hs1(k+n,:),phaseDiff,'spline');
        couplingstr = couplingstr + alpha_k*(couplingfctn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %couplingstr=couplingstr+alpha_k.*sin(theta(j)-theta(i)-(i-j)/n);
    end
    if i == m %Forcing position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        phaseDiff=theta(n+1)-theta(i);
        if(phaseDiff<0)
            phaseDiff=mod(abs(phaseDiff),0.5)*-1;
        else
            phaseDiff=mod(phaseDiff,0.5);
        end
        couplingfctn = interp1(-0.5:0.01:0.5,(H_forcing)',phaseDiff,'spline');
        couplingstr = couplingstr+alphaf*(couplingfctn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %couplingstr = couplingstr + alphaf*sin(2*pi*(theta(n+1) - theta(i)));
    end
    dy(i) = omega + couplingstr;
end
dy(n+1) = omegaf;
end

