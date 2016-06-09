function [x, l, a, G] = nn_tuned_forcing_01(pars, x0, dt, t_max)
% x = nn_01(pars, x0, dt, t_max)
%
% Lamprey segmental oscillators with smooth threshold functions
% and arbitrary intersegmental connections whose strength decreases
% exponentially as a function of connection length.
%
% Simulate for specified time using 4th-order Runge-Kutta.
%
% pars        Structure containing model parameters.
% x0          Vector of initial values--6N voltages and 1 phase of forcer.
% dt          Simulation time step.
% t_max       Simulation time.
%
% x           Variable-by-time matrix of variable values.
 
% 01-Jun-2007  Tim Kiemel
% 03-Jun-2007  Change from v01a: Take additional step when voltage crosses 0.
% 03-Jun-2007  Change from v01b: Take 2 additional steps when voltage crosses 0.
% 23-Oct-2007  Modified to include edge cell forcing
% 5 - Mar -2008  Modified to include tuning
 
%fprintf('entering nn_tuned_forcing_01\n');
%keyboard;

n       = pars.n;        % Number of oscillators.
gR      = pars.gR;       % Resting conductance for each cell type.
gT      = pars.gT;       % Tonic excitatory conductance for each cell type.
vS      = pars.vS;       % Synpatic reversal potential for each cell type.
%intra   = pars.intra;    % Intrasegmental connection parameters.
%inter   = pars.inter;    % Intersegmental connection parameters.
alpha_f = pars.alpha_f;    % Forcing strength
forcing_position = pars.forcing_position;
omegaf           = 2 * pi * pars.omegaf;
nvar=length(x0);
s = pars.s;
 
%[G l a]= G_tuned_01(pars);
G = computeIntraMatrix(pars.intra);
 
% Modified next section to include forcing
alpha_fil=alpha_f(1); 
alpha_fic=alpha_f(2); 
alpha_fel=alpha_f(3); 
alpha_fec=alpha_f(4);
 
% Edge cell connections
g = zeros(6*n,4);
g(forcing_position * 6 - 4, 1) = alpha_fel;                   % Left excitatory EC to left L.
g(forcing_position * 6 - 1, 2) = alpha_fil;                   % Left inhibitory EC to right L.
g(forcing_position * 6 - 1, 3) = alpha_fel;                   % Right excitatory EC to right L.
g(forcing_position * 6 - 4, 4) = alpha_fil;                   % Right inhibitory EC to left L.
g(forcing_position * 6 - 3, 1) = alpha_fec;                   % Left excitatory EC to left C.
g(forcing_position * 6    , 2) = alpha_fic;                   % Left inhibitory EC to right C.
g(forcing_position * 6    , 3) = alpha_fec;                   % Right excitatory EC to right C.
g(forcing_position * 6 - 3, 4) = alpha_fic;                   % Right inhibitory EC to left C.

G = [G,g];  

gR_gT = repmat(gR + gT, 2*n, 1);
gT    = repmat(gT, 2*n, 1);
vS    = repmat(vS, 2*n, 1);
vS      = [vS; 1;-1;1;-1];
 
 
dt2   = dt/2.0;
dt6   = dt/6.0;
 
n_step = round(t_max/dt);
x = zeros(nvar, n_step + 1);

x(:,1) = x0;
xc     = x0;

% for smooth h function
nv      = 6*n + 1;                    % Number of variables.
T      = 2*pi/omegaf;                 % Forcing period.
k      = round(T/dt);                 % Number of time steps per period.
r      = 2*pi/k;

for i = 2 : n_step + 1,
    x(:,i)  = rk4_step_neural(x(:,i-1), dt, dt2, dt6, gT, gR_gT, vS, G, omegaf, s);
    x(nv,i) = r*mod(i-1, k);
end