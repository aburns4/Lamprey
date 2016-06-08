function da = yourfun(t,a)
    %Define the variables
    taui=10;
    vie=1;
    el=a(8);
    ec=a(9);
    vij=[1,0,0,-1,0,0;1,0,0,-1,0,0;1,-1,0,-1,0,0;0,0,-1,0,-1,1;0,0,-1,0,0,1;0,0,-1,0,0,1];
    wij=abs(vij);
    da=zeros(9,1);
    ve=a(7);
    
    %Start calculations
    for i=1:6
        if(i==1 || i==6)
            first=ve*(vie-a(i))-(taui^-1)*a(i);
        elseif(i==2 || i==5)
            first=el*(vie-a(i))-(taui^-1)*a(i);
        else
            first=ec*(vie-a(i))-(taui^-1)*a(i);
        end
        
        second=0;
        for j=1:6
            if(a(j)>=0)
                foo=a(j);
            else
                foo=0;
            end
            second=wij(i,j)*foo*(vij(i,j)-a(i));
        end
        da(i)=first+second;
    end
end