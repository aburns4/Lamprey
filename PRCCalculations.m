clear all
close all
clc

ves=[0.005 0.0075 0.01 0.02 0.04 0.06 0.07]; %%Range of tonic drivers%% % This is eE

cycle_period = 2*pi; %%Cycle period%%

% you need to write a function that will drive the responses to use in ODE45%%
% you will find the function and the parameters in the Varkonyi paper %%

for i = 1:length(ves)
    j = length(ves);  % This is number of tonic drivers
    disp('Current level of tonic drive to E cells is')
    disp(ves(i)); disp('Please wait!');
    disp('****************************');
    ve = ves(i); % choose the ith tonic driver
    x0 = [.1 0 0 0 0 0 [ve .01 .1]];  %initial conditions for this tonic driver
    tinitial = 500; %
    numvar = 3; % number of cell types
    c1 = 10e-8; % tolerances to be used later
    c2 = 10e-8; % tolerances to be used later

    dy = 0.001;  % perturbation size

    fi{i} = 0:pi/20:cycle_period;


    %exc0=[]; % matrix of the excitation connections from column to row
    %inh0=[]; % matrix of the inhibition connections from column to row


    %------------------------length of period and initial point of period
    % the oscillator is integrated for some time until transients vanish
    [t,y]=ode45(@yourfun,[0 tinitial],x0',[]);
    %the oscillator is integrated until the beginning of a new cycle (as
    %defined in @newcycle)

    for k=1:2
        y0=y(end,:);
        [t,y]=ode45(@yourfun,[0 tinitial],y0,odeset('events',@newcycle,'Reltol',c1,'AbsTol',c2));
    end
    % length of period and initial point of a cycle are determined
    TT=t(end);
    y0=y(end,:);

    %points on one cycle
    ficycle=t*cycle_period/TT;
    ycycle=y(:,1:6);



%%%-----------------calculation of the PRC via integration of oscillators
%%%from unperturbed (y01) and perturbed (y02) initial points
    tphr=fi{i}*TT/cycle_period;
    PPRC = [];
    for index =1:numvar
        for ii=1:length(fi{i})
            if tphr(ii)==0
                t=0;
                y=y0;
            else
                [t,y]=ode45(@yourfun,[0 tphr(ii)],y0,odeset('Reltol',c1,'AbsTol',c2));
            end
            y01=y(end,:); %unperturbed initial point
            y02=y01+[zeros(1,index-1) dy zeros(1,length(x0)-index)]; %perturbed initial point
            tphr1=0;
            tphr2=0;

            for jj=1:3 %integration through 2-3 cycles
                [t1,y1]=ode45(@yourfun,[0 10*tinitial],y01,odeset('events',@newcycle,'Reltol',c1,'AbsTol',c2));
                [t2,y2]=ode45(@yourfun,[0 10*tinitial],y02,odeset('events',@newcycle,'Reltol',c1,'AbsTol',c2));
                y01=y1(end,:);
                y02=y2(end,:);
                tphr1=tphr1+t1(end);
                tphr2=tphr2+t2(end);
            end
            PPRC(ii,index)=(tphr1-tphr2)*cycle_period/TT/dy; %calculation of phase response from time differences
            percon=cycle_period/dy;
            PPRC=mod(PPRC+percon/2,percon)-percon/2; %correction for the case (tphr1-tphr2) accidentally contains a multiple of T
        end
    end

    % % %

    fi0 = ficycle;
    y0 = ycycle;

    y0=y0(:,1:6);
    fi{i}=fi{i}';

    figure(2);
    if i==1
        title ('************FIGURE 2*************');
    end
    subplot(2,j,i);

    hold on;
    plot(fi0,y0(:,1:3));
    hold off;
    subplot(2,j,j+i);
    plot(fi{i},PPRC);
    title({'FIG. 2, tonic drive =' ve});
    hold off;
    pause(0.01)
end
