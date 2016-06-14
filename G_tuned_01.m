% This function performs the tuning and outputs the H functions
% for the phase model.
% Hs1 represents the intersegmental connections.
% H_forcing represents the forcing connection.
% l represents the lengths that were successfully tuned, and
% w represents the weights given to each initial h function 
% used to determine the H functions for intersegmental and forcing
% connections.
function [l, w, Hs1, H_forcing, psi] = G_tuned_01(pars)
% pars    Structure containing model parameters.

n       = pars.n;        % Number of oscillators.
intra   = pars.intra;    % Intrasegmental connection parameters.
ra      = pars.ra;       % Ascending scaling factor.
rd      = pars.rd;       % Descending scaling factor.
la      = pars.la;       % Ascending length constant.
ld      = pars.ld;       % Descending length constant.
maxa    = pars.maxa;     % Maximum ascending length.
maxd    = pars.maxd;     % Maximum descending length.
Hfile   = pars.Hfile;    % Name of file with coupling functions.
psi_bar = pars.psi_bar;  % Desired phase lag per segment.
dw      = pars.dw;       % Weight step size to use during tuning.
%mu      = pars.mu;       % Mean of the Poisson distribution.
%seed    = pars.seed;     % Seed for Poisson distribution.
forcing_conn = pars.forcing_conn;  %Number of types of forcing connections (Set to 4 in Hprc_12_12_12)
if isfield(pars, 'g_ratio')
  g_ratio = pars.g_ratio;   % Scaling factor for intersegmental strengths.
  %fprintf('g_ratio found\n');
else
  g_ratio = 1;
end

%n_intra = size(intra, 1);    % Number of connection types.
    
load(Hfile);                 % Load file with coupling functions.
%% to test for scaling 10^-3
%H=10^(-3)*H;
% Assume that all forcing connections are given equal weight.

%GETS THE H_FORCING FROM H (WHICH CONTAINS PRCS(1-6) AND FORCING(7-10))
H_forcing = sum(H(:, end-forcing_conn + 1:end), 2);  % add all the cellular forcing connections- each column is one of the 6 cell-to-cell connections
                                                     %sum(H,2) means sum
                                                     %along the rows
%SETS H TO ONLY CONTAIN THE PRCS                                                
H=H(:,1:end-forcing_conn);

% Tune connections by adjusting weights of E-to-L and C-to-E connections.

Hsum      = sum(H, 2);  %connections summed at each phase (along rows)
%[l w Hs1] = Htune3_01(phi, 6*H(:,1), Hsum, 6*H(:,4), psi_bar, dw);
[l w Hs1 psi] = Htune3_01(phi, 6*H(:,3), Hsum, 6*H(:,2), psi_bar, dw);
% Extract weights for specified range of connection lenghts.

[m ia]    = min(abs(l + maxa));
[m id]    = min(abs(l - maxd));
l         = l(ia : id);
w         = w(ia : id, :);
Hs1       = Hs1(:, ia:id);
n_phi = size(Hs1, 1);
nl        = length(l);

% Adjust relative strength with length.
% this is for all-to-all coupling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%al = [intra * g_ratio * ra * exp((-maxa : -1)'/la); 1; intra * g_ratio * rd * exp(-(1 : maxd)'/ld)]; %coupling scheme from other models
%Hs1 = Hs1 .* repmat(al', n_phi, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%    filename = sprintf('Hfunction_data_1_23_14.mat');
%    save(filename,'H','Hs1','Hs1a','H_forcing','l','w','Hsum','Hfile','al');
%   keyboard
