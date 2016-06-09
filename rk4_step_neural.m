% Fourth-order Runge-Kutta for the neural model
function x = rk4_step_neural(xp, dt, dt2, dt6, gT, gR_gT, vS, G, wf, s)

nv = length(xp);    % Number of variables.
k  = 1 : nv - 1;    % Indices of CPG variables.
f  = zeros(nv, 4);

for j = 1 : 4,
  if j == 1
    x = xp;                         % wf
  elseif j == 2
    x = xp + 0.5*dt*f(:,1);
  elseif j == 3
    x = xp + 0.5*dt*f(:,2);
  else
    x = xp + dt*f(:,3);
  end

  v1      = x(k);            % Voltages of CPG cells.
  v2      = sin(x(nv));      % Voltage of left edge cells.

  v       = [v1; v2; v2; -v2; -v2];

  %smooth
  hv = s * log(1 + exp(v / s));
  
  
  f(k,j)  = gT + G*(hv.*vS) - (gR_gT + G*hv).*v1;

  f(nv,j) = wf;
end

x = xp + (dt6)*(f(:,1) + 2*(f(:,2) + f(:,3)) + f(:,4)); 