a = 1;
x = zeros(1,15);
y = zeros(1,15);
z = zeros(1,15);
for k=1:15
    h = 10.^-k;
    x(k) = h;
    g1 = (exp(a+h)-exp(a))/h;
    g2 = (exp(a+h)-exp(a-h))/(2*h);
    y(k) = abs(exp(a)-g1)/exp(a);
    z(k) = abs(exp(a)-g2)/exp(a);
end
figure(1)
plot(x,y,'r',x,z,'b')
figure(2)
loglog(x,y,'r',x,z,'b')


