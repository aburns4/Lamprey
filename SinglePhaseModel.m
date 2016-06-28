omega = 1;
dc = 0.36;
Aa = 1.0;
Ad = 10.0;
lambda_a = 40;
lambda_d = 5;
n = 280;
lag = 1/n;

alpha_f = 4;
omega_f = 1;
fPos = 0;

theta_0 = [0:-lag:-1+lag 0];

t_0 = 0;
t_f = 10;
dt = 0.001;
N = floor((t_f-t_0)/dt);

t = (t_0:dt:t_f);
theta = zeros(N,n+1);
theta(1,:) = theta_0;
theta_dot = zeros(1,n+1);

for timeStep = 1:N
    for i = 1:n
        coupling = 0;
        for j = 1:n
            if(i-j < 0)
                s = Aa*exp(-abs(i-j)/lambda_a);
            elseif(i-j > 0)
                s = Ad*exp(-abs(i-j)/lambda_d);
            else
                s = 0;
            end
            phase_lag = (i-j)/n;
            coupling = coupling + s*sin(theta(timeStep,j)-theta(timeStep,i)-phase_lag);
        end
        
        forcing = 0;
        if(i == fPos)
            forcing = alpha_f*sin(theta(timeStep,n+1)-theta(timeStep,i));
        end
        
        theta_dot(i) = omega + coupling + forcing;
    end
    theta_dot(n+1) = omega_f;
    theta(timeStep+1,:) = theta(timeStep,:) + theta_dot*dt;
end