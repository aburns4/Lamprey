omega = 1;
dc = 0.36;
Aa = 1.0;
Ad = 10.0;
lambda_a = 40;
lambda_d = 5;
n = 280;
lag = 0.01;

alpha_f = 4;
omega_f = 1.0;
fPos = 150;
fSide = 1;

theta_0 = [0:-lag:-1+lag 0];

t_0 = 0;
t_f = 4;
dt = 0.001;
N = floor((t_f-t_0)/dt);

t = (t_0:dt:t_f);
theta = cell(1,2);
theta{1} = zeros(N,n+1);
theta{2} = 0.5*ones(N,n+1);
%theta{1}(1,:) = theta_0;
%theta{2}(1,:) = theta_0;
theta_dot = zeros(1,n+1);

side_strength = 10*Ad*exp(-1/lambda_d);

for timeStep = 1:N
    for side = 1:2
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
                phase_lag = (i-j)*lag;
                coupling = coupling + s*sin(2*pi*(theta{side}(timeStep,j) - theta{side}(timeStep,i) - phase_lag));
            end
            
            forcing = 0;
            if(side == fSide)
                if(i == fPos)
                    forcing = alpha_f*sin(theta{side}(timeStep,n+1)-theta{side}(timeStep,i));
                end
            end
            
            if(side == 1)
                cross_conn = side_strength*sin(2*pi*(theta{2}(timeStep,i)-theta{1}(timeStep,i))+pi);
            else
                cross_conn = side_strength*sin(2*pi*(theta{1}(timeStep,i)-theta{2}(timeStep,i))+pi);
            end
                
            theta_dot(i) = omega + coupling + forcing + cross_conn;
        end
        theta_dot(n+1) = omega_f;
        theta{side}(timeStep+1,:) = theta{side}(timeStep,:) + theta_dot*dt;
    end
end

% segments = 1:n;
% close all;
% figure(1);
% for timeStep = 1:N
%     if(mod(timeStep,10)==0)
%         side1 = sin(2*pi*theta{1}(timeStep,segments));
%         side2 = sin(2*pi*theta{2}(timeStep,segments));
%         
%         plot(segments,side1,segments,side2);
%         axis([0 n+1 -1 1]);
%         title(sprintf('Time is %0.3f',timeStep*dt));
%         pause(0.01);
%     end
% end

activation = cell(1,2);
activation{1} = zeros(N,n);
activation{2} = zeros(N,n);

for side = 1:2
    for timeStep = 1:N
        for i = 1:n
            if(sin(2*pi*theta{side}(timeStep,i)) >= sin(2*pi*omega*(0.25*omega+dc/2)))
                activation{side}(timeStep,i) = 1;
            else
                activation{side}(timeStep,i) = 0.1;
            end
        end
    end
end     
    
segments = 1:n;
close all;
figure(1);
for timeStep = 1:N
    if(mod(timeStep,10)==0)
        side1 = activation{1}(timeStep,segments);
        side2 = -activation{2}(timeStep,segments);
        
        plot(segments,side1,segments,side2);
        axis([0 n+1 -1.1 1.1]);
        title(sprintf('Time is %0.3f',timeStep*dt));
        pause(0.01);
    end
end
    
    
    
    