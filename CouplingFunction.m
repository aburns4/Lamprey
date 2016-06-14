%% GET PRC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 1;                               % Number of oscillators.
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
%%%%%
ra      = 0.002;         % Amplitude of ascending strength function.
rd      = 0.0005;        % Amplitude of descending strength function.
g_ratio = 20;            % Scaling factor for intersegmental strenghts.
la      =  3;            % Length constant of ascending strength function.
ld      = 10;            % Length constant of descending strength function.
maxa    = n-1;            % Maximum ascending length. All to all.
maxd    = n-1;            % Maximum descending length. All to all.     
%%%%%
Sym = 1;
     
mu      = [];            % Average number of connections per segment
seed    = 0;             % Seed for generating random connections.
%%%%%
Hfile   = 'H_prc_11314';%'H_forcing06_001';     % Name of file with coupling functions.
%%%%%
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
                      'omegaf',omegaf,'forcing_conn',forcing_conn,'Sym',Sym,...%);
                        'ra',ra,'rd',rd,'la',la,'ld',ld,'g_ratio',g_ratio,...
                            'maxa',maxa,'maxd',maxd, 'Hfile', Hfile);%, ...
                 

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GET H_C FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hprc
% This function takes as input the phase shifts resulting from perturbing
% the voltages of each cell in an oscillator, at each phase. It produces
% functions representing the strength of intersegmental ascending and
% descending connections of varying lengths. For instance, for a chain of 3
% oscillators, the output Hprc (describing ascending connections)
% would consist of two rows - the first would represent ascending 
% connections of length 1, and the second would represent ascending 
% connections fo length 2.
% The descending connections are formatted similarly as H_desc.

phi = -0.5:0.01:0.5;
n_phi = length(phi);
n_phi2 = (n_phi + 1)/2;

intra = pars.intra;
intra(:,3) = 2 * intra(:,3); % Double the connection strengths because there are two of each type.

s = pars.s;

n_cell = 6 * pars.n;
pars.n = 3;               % use a chain of 3 neural oscillators
vS = [pars.vS; pars.vS];  %synaptic reversal potentials (1 or -1)

N = size(intra,1);

x = x';  % voltage of each cell during last oscillation
P = P'; %matrix of phase shifts (PRCs)
x = x(:,1:100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the H function for ascending connections. 
% Each row of asc has the form (from, to, strength).
H_inter = zeros(n_phi,N);
connLength = 1;

for connNum = 1:N
    fromCell = intra(connNum,1); % these are like i,j,g in neural model
    toCell = intra(connNum,2);
    strength = intra(connNum,3);
    for delphi = 0:n_phi-1                  % number of phases n_phi
        points = zeros(n_phi-1,1);
        for fromIndex = 1:n_phi-1             % the sender
            toIndex = mod(fromIndex - delphi - 1, n_phi-1) + 1;
            sender = x(1:6, fromIndex);   %voltage of each of six cells in the sender oscillator
            receiver = x(1:6, toIndex);   %voltage of cells for receiving oscillator
            
            hval = s * log(1 + exp(sender(fromCell,1) / s));        %smooth h function
            %hval = max(sender(fromCell,1), 0);                    
            vSval = vS(fromCell, 1);
            change = strength * hval * (vSval - receiver(toCell, 1));  %strength of connection*smoothing function* driving force (difference between reveral potential and cell voltage)
            %this is the connection between the two cells in the neural model                                                        
            points(toIndex, 1) = P(toCell,toIndex) * change;
        end
        %H_inter(delphi + 1, connNum) = (1 / (n_phi - 1)) * (1/2) * dot(points, [1; repmat(2, n_phi - 2, 1); 1]);
        H_inter(delphi + 1, connNum) = mean(points);
        % We use the trapezoidal approximation of the integral. The line above
        % averages the effects of a specific intersegmental connection
        % for a specific phase shift delphi, over all phases.
    end
end
H_inter = [H_inter(n_phi2:end,:); H_inter(2:n_phi2,:)]; %rearrange so phases go from -0.5 to 0.5
%H_inter = H_inter * T;       % Why T?

titleList = ['E to L'; 'E to C';
             'L to C'; 'C to E';
             'C to L'; 'C to C'];

phi1 = -0.5:0.01:0.5;

figure(1)
for i = 1 : 3,
  for j = 1 : 2,
    k = 2*(i - 1) + j;
    subplot(3,2,k);
    plot(phi1, H_inter(:,k), '-', 'LineWidth', 2);
    title(titleList((i-1)*2 + j,:), 'FontSize', 18);
    xlabel('Relative phase', 'FontSize', 12);
    ylabel('Connection strength', 'FontSize', 8);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the H function for forcing.
H_forcing = zeros(n_phi, 4);
fraction = 1;
alpha_fil =15 * fraction;
alpha_fic =15 * fraction;
alpha_fel =15 * fraction;
alpha_fec =15 * fraction;
alpha_f = [alpha_fil, alpha_fic, alpha_fel, alpha_fec];
pars.alpha_f = 2 * alpha_f;     % double the connection strength because there are two of each type of connection.
%pars.g_ratio = 20;
% The edge cell connection excites the same sided cells (so Vsyn_ec =
% 1 for thes cells) and inhibits the opposite sided cells (so Vsyn_ec = -1);
Vsyn_ec = [-1, -1,  1,  1];
jIndices = [5, 6, 2, 3];

% init is the forcing phase associated with the oscillator phase = 0.
init = 0;
forcerPhase = init: (2*pi/100): init + 2*pi;
for j = 1:4
    for delphi = 0:n_phi-1 
        points = zeros(n_phi-1,1);
        for fromIndex = 1:n_phi-1
            toIndex = mod(fromIndex - delphi - 1, n_phi - 1) + 1;
            G = alpha_f(1, j); % double the connection strength because there are two of each type of connection.
            v_ec = (-1)^ ceil((j) / 2) * sin(forcerPhase(fromIndex));
            Vsyn = Vsyn_ec(1,j);
            vij = x(jIndices(j), toIndex);
            change = G * (s * log(1 + exp(v_ec / s))) * (Vsyn - vij);       %smooth h
            points(toIndex, 1) = P(toCell,toIndex) * change;
        end
        %H_forcing(delphi + 1, j) = (1 / (n_phi - 1)) * (1/2) * dot(points, [1; repmat(2, n_phi - 2, 1); 1]);
        H_forcing(delphi + 1, j) = mean(points);
    end
end

H_forcing = [H_forcing(n_phi2:end,:); H_forcing(2:n_phi2,:)];
%H_forcing = H_forcing * T;

titleList = ['EC to L inhibitory'; 'EC to C inhibitory';
             'EC to L excitatory'; 'EC to C excitatory'];
figure(2);
for i = 1 : 2,
  for j = 1 : 2,
    k = 2*(i - 1) + j;
    subplot(2,2,k);
    plot(phi1, H_forcing(:,k), '-', 'LineWidth', 2);
    grid on;
    title(titleList(2*(i-1) + j, :), 'FontSize', 18);
    xlabel('Relative phase', 'FontSize', 12);
    ylabel('Connection strength', 'FontSize', 12);
  end
end

H = [H_inter, H_forcing];

%Our H_forcing_final???
H_forcing_final=2*sum(H_forcing,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GET W(weights of connections)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set Number of oscillators needed
pars.n       = 5;                  % Number of oscillators
n = pars.n;
pars.maxa    = pars.n-1;            % Maximum ascending length. All to all.
pars.maxd    = pars.n-1;            % Maximum descending length. All to all. 


%Run G_tuned_01 to get w(weights of connections) and l(distance)
%Hs1 -> H_c, H_forcing -> H_f (forcing function)
%w columns are L to C, other, E to C
[l, w, Hs1, H_forcing, psi] = G_tuned_01(pars);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GET CORRECT VALUES FOR WEIGHTS AND CALCULATE H_R (COUPLING) FUNCTION

numRows = size(w,1);

for i=1:numRows
    if(w(i,1)==0)
        break;
    end
end

other=w(:,2);

LtoC = zeros(numRows,1);
LtoC(1:i-1) = 6*w(1:i-1,1)+other(1:i-1);
LtoC(i:numRows) = other(i:numRows);

EtoC = zeros(numRows,1);
EtoC(1:i-1) = other(1:i-1);
EtoC(i:numRows) = 6*(w(i:numRows,3)) + other(i:numRows);

% figure(3);
% scatter(-(n-1):(n-1),other);
% hold on;
% scatter(-(n-1):(n-1),EtoC);
% scatter(-(n-1):(n-1),LtoC);
% hold off;

%matrix of values [EtoL, EtoC, LtoC, CtoE, CtoL, CtoC]
connTypes=[other,EtoC,LtoC,other,other,other];
H_inter=H_inter';

%Multiply H_inter by connTypes to compute H_r
%This is (4) is Massarelli, 2016
H_r=connTypes*H_inter;

save('CouplingFunction.mat','H_r','H_forcing');