function [theta,STRENGTH] = allToall(t,theta,n,w,forcingTheta,forcingOmega,alphaf)
    STRENGTH=zeros(n,n);
    forcing=zeros(n,1);
    for i=1:n
        sum=0;
        for j=1:n-1
            k=mod(i+j-1,n)+1;
            sum=sum+fstrength(i,k,n,theta(k),theta(i));
            STRENGTH(i,k)=fstrength(i,k,n,theta(k),theta(i));
        end
        if(forcingOmega(i)~=0)
            forcing(i)=alphaf*sin(forcingTheta(i)-theta(i));
        end
        theta(i)=w+sum;
    end
    theta=theta+forcing;
end
