function [tb, xb] = nn_tuned_bursts_01(pars, x0, dt, ib, nb_min, t_max)
% [tb, xb] = nn_tuned_bursts_01(pars, x0, dt, ib, nb_min, t_max)
%
% Lamprey segmental oscillators with smooth threshold functions
% and arbitrary intersegmental connections whose strength decreases
% exponentially as a function of connection length.
%
% Compute specified number of zero crossings for the first phase model
% oscillator using 4th-order Runge-Kutta.
%
% pars        Structure containing model parameters.
% x0          Vector of initial values.
% dt          Simulation time step.
% ib          Index of reference cell.
% nb_min      Minimum number of bursts times per segment to compute.
% t_max       Return early if simulation time reaches t_max.
%
% tb          Cell array of burst times.
% xb          Variable values at last burst time.

% 03-Jun-2007  Tim Kiemel.
% 30-Jul-2007  Change from v1: Added ib input.
% 17-Jul-2012  Added tuning

n       = pars.n;        % Number of oscillators.
gR      = pars.gR;       % Resting conductance for each cell type.
gT      = pars.gT;       % Tonic excitatory conductance for each cell type.
vS      = pars.vS;       % Synpatic reversal potential for each cell type.
%%%%%
%intra   = pars.intra;    % Intrasegmental connection parameters.
%inter   = pars.inter;    % Intersegmental connection parameters.
%%%%%
alpha_f = pars.alpha_f;    % Forcing strength
forcing_position = pars.forcing_position;
omegaf           = 2 * pi * pars.omegaf;
nvar=length(x0);
s = pars.s;
%%%%%
%[G l a]= G_tuned_01(pars);
%%%%%
G = computeIntraMatrix(pars.intra);
%%%%%

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
%keyboard

gR_gT = repmat(gR + gT, 2*n, 1);
gT    = repmat(gT, 2*n, 1);
vS    = repmat(vS, 2*n, 1);
vS      = [vS; 1;-1;1;-1];
 
dt2   = dt/2.0;
dt6   = dt/6.0;
 
%keyboard;
n_step = round(t_max/dt);
xc = x0;

nb     = zeros(n, 1);
nb_max = (nb_min + 100)*ones(n, 1);  % Allocate memory for nb_max burst times.
tb     = cell(n, 1);
for i = 1 : n,
    tb{i} = zeros(nb_max(i), 1);
end

for i = 1 : n_step,
    xp  = xc;
    xc  = rk4_step_neural(xp, dt, dt2, dt6, gT, gR_gT, vS, G, omegaf, s);
    
    j   = find(xp < 0 & xc >= 0 | xp > 0 & xc <= 0);
    nj  = size(j, 1);
    if nj > 0
        dtc     = dt * (xp(j) ./ (xp(j) - xc(j)));
        [dtc k] = sort(dtc);
        j       = j(k);
        
        dt_sum  = 0;
        for k = 1 : nj,
            jk     = j(k);
            
            dtc1   = dtc(k) - dt_sum;
            x1     = rk4_step_neural(xp, dtc1, dtc1/2, dtc1/6, gT, gR_gT, vS, G, omegaf, s);
            
            dtc2   = dtc1 * xp(jk) / (xp(jk) - x1(jk));
            x2     = rk4_step_neural(xp, dtc2, dtc2/2, dtc2/6, gT, gR_gT, vS, G, omegaf, s);
            
            if mod(jk, 6) == ib & xc(jk) >= 0         % Burst time for left E.
                dtc3   = dtc2 * xp(jk) / (xp(jk) - x2(jk));
                xp     = rk4_step_neural(xp, dtc3, dtc3/2, dtc3/6, gT, gR_gT, vS, G, omegaf, s);
                
                dt_sum = dt_sum + dtc3;
                
                p       = floor((jk + 5)/6);           % Segment index.
                nb(p)   = nb(p) + 1;
                
                if nb(p) > nb_max(p)                   % If necessary,
                    tb{p}     = [tb{p}; zeros(100,1)];   % allocate more memory for
                    nb_max(p) = nb_max(p) + 100;         % burst times.
                end
                
                tb{p}(nb(p)) = (i - 1)*dt + dt_sum;
                
                %if min(nb) == nb_min
                if nb(1) == nb_min
                    xb = xp;
                    for r = 1 : n,
                        tb{r} = tb{r}(1 : nb(r));
                    end
                    return;
                end
            else
                xp     = x2;
                dt_sum = dt_sum + dtc2;
            end
        end
        
        dtc3 = dt - dt_sum;
        xc   = rk4_step_neural(xp, dtc3, dtc3/2, dtc3/6, gT, gR_gT, vS, G, omegaf, s);
    end
end

xb = xc;
for r = 1 : n,
    tb{r} = tb{r}(1 : nb(r));
end