function da = networkFunc(t,a,tDrives,v_ie,T_i,v_ij,w_ij,n)

    da = zeros(6*n,1);
    
    for i = 1:n
        for j = 1:6
            firstTerm = 0;
            if(j==1 || j==4) %E cells
                firstTerm = tDrives(1)*(v_ie-a(6*(i-1)+j))-(1/T_i)*a(6*(i-1)+j);
            elseif(j==2 || j==5) %L cells
                firstTerm = tDrives(2)*(v_ie-a(6*(i-1)+j))-(1/T_i)*a(6*(i-1)+j);
            elseif(j==3 || j==6)%C cells
                firstTerm = tDrives(3)*(v_ie-a(6*(i-1)+j))-(1/T_i)*a(6*(i-1)+j);
            end
            if(n == 1)
                lowerBound = 1;
                upperBound = 1;
            elseif(i-1 < 1) %If "i" is on the beginning edge of the chain
                lowerBound = 1;
                upperBound = i+1;
            elseif(i+1 > n) %If "i" is on the 
                lowerBound = i-1;
                upperBound = n;
            else
                lowerBound = i-1;
                upperBound = i+1;
            end
            
            secondTerm = 0;
            for k = lowerBound:upperBound
                for l = 1:6
                    if(i-k == 0)
                        w_ik = w_ij(1);
                    elseif(i-k == -1)
                        w_ik = w_ij(2);
                    elseif(i-k == 1)
                        w_ik = w_ij(3);
                    end
                    
                    if((l==1 && j==2) || (l==4 && j==5) || (l==1 && j==3) || (l==4 && j==6)) %E to L and E to C excitatory connections
                        if(a(6*(k-1)+l)<0)
                            f_a=0;
                        else
                            f_a=a(6*(k-1)+l);
                        end
                        secondTerm = secondTerm + w_ik*f_a*(v_ij(1)-a(6*(i-1)+j));
                    elseif((l==2 && j==3) || (l==5 && j==6) || (l==3 && (j==4 || j==5 || j==6)) || (l==6 && (j==1 || j==2 || j==3))) %L to C and C to E,L,C inhibitory connections
                        if(a(6*(k-1)+l)<0)
                            f_a=0;
                        else
                            f_a=a(6*(k-1)+l);
                        end
                        secondTerm = secondTerm + w_ik*f_a*(v_ij(2)-a(6*(i-1)+j));
                    end
                end
            end
            da(6*(i-1)+j) = firstTerm + secondTerm;
        end
    end
end