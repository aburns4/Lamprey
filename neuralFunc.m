function dv = neuralFunc(t,v,n,m,G_R,G_T,G_0,V_syn,G_f,V_synec,sigma,alpha_f,alpha_r,omega_f)
    
    dv=zeros((6*n)+1,1);

    for i = 1:n
        result = 0;
        for j = 1:6
            result = -G_R*v(6*(i-1)+j);
            if(j == 1 || j == 4) %E cell
                result = result + G_T(1)*(1-v(6*(i-1)+j));
            elseif(j == 2 || j == 5) %L cell
                result = result + G_T(2)*(1-v(6*(i-1)+j));
            else %C cell
                result = result + G_T(3)*(1-v(6*(i-1)+j));
            end
            
            %Double Sum Part of Equation
            doublesum = 0;
            for k = 1:n
                for l = 1:6
                    
                    %Setting up index for alpha_r (strength value matrix)
                    r=i-k;
                    r=r+n;
                    
                    alpha_r(n,:)=ones(1,3);
                    
                    if((l==2 && j==3) || (l==5 && j==6)) %L to C connection
                        strength=alpha_r(r,1);
                        doublesum = doublesum + strength*G_0(1)*(sigma*log(1+exp(v(6*(k-1)+l)/sigma)))*(V_syn(2)-v(6*(i-1)+j));
                    elseif((l==1 && j==3) || (l==4 && j==6)) %E to C connection
                        strength=alpha_r(r,2);
                        doublesum = doublesum + strength*G_0(2)*(sigma*log(1+exp(v(6*(k-1)+l)/sigma)))*(V_syn(1)-v(6*(i-1)+j));
                    elseif((l==1 && j==2) || (l==4 && j==5)) %E to L connection
                        strength=alpha_r(r,3);
                        doublesum = doublesum + strength*G_0(2)*(sigma*log(1+exp(v(6*(k-1)+l)/sigma)))*(V_syn(1)-v(6*(i-1)+j));
                    else %All C to other connection
                        strength=alpha_r(r,3);
                        doublesum = doublesum + strength*G_0(2)*(sigma*log(1+exp(v(6*(k-1)+l)/sigma)))*(V_syn(2)-v(6*(i-1)+j));
                    end
                end
            end
            result = result+doublesum;
            
            %Forcing Term
            if(i == m)
                forcingSum = 0;
                for s = 1:2 %For each EC side
                    forcingSum = G_f*(sigma*log(1+exp(((-1)^s*sin(2*pi*v(end)))/sigma)));
                    if(s==1 && (j==2 || j==3)) %Left side, EC to L or C
                        forcingSum = forcingSum*(V_synec(1)-v(6*(i-1)+j));
                    elseif(s==1 && (j==5 || j==6)) %Left side, EC to oppposite L or C
                        forcingSum = forcingSum*(V_synec(2)-v(6*(i-1)+j));
                    elseif(s==2 && (j==5 || j==6)) %Right side, EC to L or C
                        forcingSum = forcingSum*(V_synec(1)-v(6*(i-1)+j));
                    elseif(s==2 && (j==2 || j==3)) %Right side, EC to opposite L or C
                        forcingSum = forcingSum*(V_synec(2)-v(6*(i-1)+j));
                    else
                        forcingSum = 0;
                    end
                end
                result = result + alpha_f*forcingSum;
            end
            
            dv(6*(i-1)+j)=result;
        end
    end
    dv(end)=omega_f;
end