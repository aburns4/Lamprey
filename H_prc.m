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

load prc_121212.mat;

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

save('H_prc_11314.mat', 'H', 'H_inter', 'H_forcing', 'pars', 'P', 'phi', 'n_phi', 'x', 'H_forcing_final');

