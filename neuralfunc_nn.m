function [dv] = neuralfunc_nn(v_0,n,m,G_R,G_T,G_0,V_syn,G_f,V_synec,sigma,alpha_f,alpha_r,omega_f)

dv=zeros(6*n+1,1);

for i = 1:n
    for j=1:6
        if n==1
            lowerbound=1;
            upperbound=1;
        elseif i==1 %if first segment
            lowerbound=1;
            upperbound=2;
            
        elseif i==n %if last segemnt
            lowerbound=n-1;
            upperbound=n;
            
        else %all other segments
            lowerbound=i-1;
            upperbound=i+1;
        end
        
        result=0;
        result = -G_R*v_0(6*(i-1)+j);
        if(j == 1 || j == 4) %E cell
            result = result + G_T(1)*(1-v_0(6*(i-1)+j));
        elseif(j == 2 || j == 5) %L cell
            result = result + G_T(2)*(1-v_0(6*(i-1)+j));
        else %C cell
            result = result + G_T(3)*(1-v_0(6*(i-1)+j));
        end
        
        %Double Sum Part of Equation
        doublesum = 0;
        for k = lowerbound:upperbound
            for l = 1:6
                
                %Setting up index for alpha_r (strength value matrix)
                r=i-k;
                r=r+2;
                strength=alpha_r(r);
                
                
                if((l==2 && j==3) || (l==5 && j==6)) %L to C connection
                    doublesum = doublesum + strength*G_0(1)*(sigma*log(1+exp(v_0(6*(k-1)+l)/sigma)))*(V_syn(2)-v_0(6*(i-1)+j));
                elseif((l==1 && j==3) || (l==4 && j==6)) %E to C connection
                    doublesum = doublesum + strength*G_0(2)*(sigma*log(1+exp(v_0(6*(k-1)+l)/sigma)))*(V_syn(1)-v_0(6*(i-1)+j));
                elseif((l==1 && j==2) || (l==4 && j==5)) %E to L connection
                    doublesum = doublesum + strength*G_0(2)*(sigma*log(1+exp(v_0(6*(k-1)+l)/sigma)))*(V_syn(1)-v_0(6*(i-1)+j));
                elseif((l==3 && j==4) || (l==6 && j==1)) %C to E connection
                    doublesum = doublesum + strength*G_0(2)*(sigma*log(1+exp(v_0(6*(k-1)+l)/sigma)))*(V_syn(2)-v_0(6*(i-1)+j));
                elseif((l==3 && j==5) || (l==6 && j==2)) %C to L connection
                    doublesum = doublesum + strength*G_0(2)*(sigma*log(1+exp(v_0(6*(k-1)+l)/sigma)))*(V_syn(2)-v_0(6*(i-1)+j));
                elseif((l==3 && j==6) || (l==6 && j==3)) %C to C connection
                    doublesum = doublesum + strength*G_0(2)*(sigma*log(1+exp(v_0(6*(k-1)+l)/sigma)))*(V_syn(2)-v_0(6*(i-1)+j));
                end
            end
        end
        result = result+doublesum;
        
        %Forcing Term
        if(i == m)
            finalSum = 0;
            for s = 1:2 %For each EC side
                tempSum = 0;
                forcingSum = G_f*(sigma*log(1+exp(((-1)^s*sin(2*pi*v_0(end)))/sigma)));
                if(s==1 && (j==2 || j==3)) %Left side, EC to L or C
                    tempSum = tempSum + forcingSum*(V_synec(1)-v_0(6*(i-1)+j));
                elseif(s==1 && (j==5 || j==6)) %Left side, EC to oppposite L or C
                    tempSum = tempSum + forcingSum*(V_synec(2)-v_0(6*(i-1)+j));
                elseif(s==2 && (j==5 || j==6)) %Right side, EC to L or C
                    tempSum = tempSum + forcingSum*(V_synec(1)-v_0(6*(i-1)+j));
                elseif(s==2 && (j==2 || j==3)) %Right side, EC to opposite L or C
                    tempSum = tempSum + forcingSum*(V_synec(2)-v_0(6*(i-1)+j));
                else
                    tempSum = 0;
                end
                finalSum = finalSum + tempSum;
            end
            result = result + alpha_f*finalSum;
        end
        dv(6*(i-1)+j)=result;
    end
end

dv(end)=omega_f;
end