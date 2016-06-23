%First open the downloaded folder, then run

A = dlmread('fort.55554');
B = reshape(A(:,2),319,20);
 
for i=1:size(B,1)
    plot(sin(2*pi*B(i,:)))
    axis([1 size(B,2) -1 1])
    pause(.01)
end