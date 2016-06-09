% Computes the PRCs (phase response curves) for the tuned neural model.
function [omega, prc] = get_prc()

n = 1;                                % Number of oscillators.
forcing_position = 1;                 % There is no forcing in this run.
n_cell = 6;                           % Number of cells per oscillators.

gR      = [ 3.500;  3.500;  3.500];   % Resting conductance for each cell type.
gT      = [ 0.875;  0.350;  3.500];   % Tonic excitatory conductance for
                                      % each cell type.
vS      = [ 1    ; -1    ; -1    ];   % Synpatic reversal potential for
                                      % each cell type.
s       = 0.05;                       % Scale for smoothing threshold function.
fprintf('Smooth h, s = %d\n', s);

intra = [ 1  2  35.0;
          1  3  35.0;
          2  3  15.0;   % change to 15
          3  4  35.0;
          3  5  35.0;
          3  6  35.0];

%ra      = 0.002;         % Amplitude of ascending strength function.
%rd      = 0.0005;        % Amplitude of descending strength function.
%g_ratio = 20;            % Scaling factor for intersegmental strenghts.
%la      =  3;            % Length constant of ascending strength function.
%ld      = 10;            % Length constant of descending strength function.
%maxa    = n-1;            % Maximum ascending length. All to all.
%maxd    = n-1;            % Maximum descending length. All to all.     

Sym = 1;
     
mu      = [];            % Average number of connections per segment
seed    = 0;             % Seed for generating random connections.

%Hfile   = 'H_forcing06_001';     % Name of file with coupling functions.
psi_bar = 0.010;                 % Desired phase lag per segment.
Tuning  = 1;                     % Tuning method.
dw      = 0.001;                 % Weight step size to use during tuning.
forcing_conn     = 4;

omegaf = 1;                      % The forcing frequency.
alpha_f = [0;0;0;0];             % No edge cell connections. 

pars = struct('n',n,'gR',gR,'gT',gT,'vS',vS,'s',s, ...
              'intra',intra,'alpha_f',alpha_f, ...
                'psi_bar',psi_bar, ...
                  'Tuning',Tuning,'dw',dw','mu',mu,'seed',seed, ...
                    'forcing_position', forcing_position, ...
                      'omegaf',omegaf,'forcing_conn',forcing_conn,'Sym',Sym);
                  
               %'ra',ra,'rd',rd,'la',la,'ld',ld,'g_ratio',g_ratio,...
               %'maxa',maxa,'maxd',maxd. 'Hfile', Hfile, ...
                 

% Consider all phases from 0 to 2*pi, using increments of size 2*pi / 100.
% We scale the phases from 0 to 2*pi down to 0 to 1.
dphi = 0.01;
phi = (0 : dphi : 1);
n_phi = length(phi);
n_phi2 = (n_phi + 1)/2;

dx = 1e-5;              % The size of our perturbation.

x0      = [ones(3,n); -ones(3,n)];
x0      = x0(:);
x0      = [x0(:);0];

dtn   = 0.02;      % Simulation time step
t_max = 100;            % Simulation time.

% Simulate without forcing to determine T and the set of voltages over one
% period.
ib = 1;
nb1 = 10;
nb2 = 4;
[t x0] = nn_tuned_bursts_01(pars, x0, dtn, ib, nb1, t_max);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Find the period T and frequency omegaf of the chain without 
% any forcing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the negative-to-positive voltage crossings of the left E cell.

T = t{1}(end) - t{1}(end-1);
pars.T = T;
pars.omega0 = 1 / T;        % the natural frequency
stepsPerPeriod = round(T / dtn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: For each phase, determine the voltage of each cell, during the
% last full oscillation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt     = dphi*T/ceil(dphi*T/dtn);  % dphi*T is multiple of time step.
k      = round(dphi*T/dt);

x = nn_tuned_forcing_01(pars,x0,dt,T);
x      = x(:,1:k:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: For each phase, perturb each of the oscillator's voltages by dx, 
% one at a time. Then simulate a for time nT, where T is the period of
% oscillation. Determine the resulting phase shift, and save that in a
% matrix called prc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = zeros(n_phi, n_cell/2);

total_steps = t_max / dtn;
numPeriods = floor(total_steps / stepsPerPeriod);   
% numPeriods = number of periods being simulated when the inputs to
%              nn_tuned_forcing_01 are dtn and t_max.

for i = 1:n_phi
    for j = 1:3
        x0 = x(:,i);
        % Perturb one of the voltages by dx.
        x0(j) = x0(j) + dx;
        
        % Then simulate a fixed number of periods to
        % determine the phase shift resulting from
        % the perturbation.
        t = nn_tuned_bursts_01(pars, x0, dt, ib, nb2, t_max);
        
        P(i,j) = -t{1}(end)/T - phi(i);      
        fprintf('prc(%d, %d) = %20.14f\n', j, i, (P(i,j) - round(P(i,j))) / dx);
    end
end

P = P - round(P);     % All prc values are now between -0.5 and 0.5.
P = P / dx;
P      = [P  [P(n_phi2:end,:); P(2:n_phi2,:)] ];

x = x';

save('prc_test.mat', 'P', 'T', 'pars', 'x', 'n_phi', 'phi', 'dx');