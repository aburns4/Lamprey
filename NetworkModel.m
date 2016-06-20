tDrives = [0.01, 0.01, 0.1];
T_i = 10;
v_ie = 1;
v_ij = [1, -1];
w_ij = [1, 0.05, 0.02];
dt = 0.001;
n=5;
a_0 = zeros(6*n,1); %[-0.5, -0.5, 0, 0 ,0.3, 0.5]; 
options=odeset('RelTol',1e-4,'AbsTol',1e-7);
[T, a] = ode45(@networkFunc, (0:dt:100), a_0, options, tDrives, v_ie, T_i, v_ij, w_ij, n);
plot(T, a);