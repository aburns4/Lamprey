clear all
close all
clc

ves=[.005 .0075 .01 .02 .04 .06 .07]; % This is eE

for i = 1:length(ves)
    
    
    j = length(ves);  % This is number of tonic drivers
    
    ve = ves(i); % choose the ith tonic driver
    disp('Current level of tonic drive to E cells is')
    disp(ves(i)); disp('Please wait!');
    disp('****************************');
    x0 = [.1 0 0 0 0 0 [ve .01 .1]];  %initial conditions for this tonic driver
    tinitial = 500; % 
    numvar = 3; % number of cell types
    c1 = 10e-8; % tolerances to be used later
    c2 = 10e-8; % tolerances to be used later
    
    dy = 0.001;  % step size
    
    fi{i} = 0:pi/20:2*pi; % this breaks up the period into intervals, maybe used for time stepping? 
    
   
    exc0=[0 1 1 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 1 1 0]'; % matrix of the excitation connections from column to row
    inh0=[0 0 0 0 0 0;0 0 1 0 0 0;0 0 0 1 1 1;1 1 1 0 0 0;0 0 0 1 0 0;0 0 0 0 0 0]'; % matric of the inhibition connections from column to row
   
    
    %------------------------length of period and initial point of period
    % the oscillator is integrated for some time until transients vanish
    [t,y]=ode45(@mywillfun,[0 tinitial],x0');
    %the oscillator is integrated until the beginning of a new cycle (as
    %defined in @newcycle)
       
    for k=1:2
        y0=y(end,:);
        [t,y]=ode45(@mywillfun,[0 tinitial],y0,odeset('events',@newcycle,'Reltol',c1,'AbsTol',c2));
    end
    % length of period and initial point of a cycle are determined
    TT=t(end);
    y0=y(end,:);
    
    %points on one cycle
    ficycle=t*2*pi/TT;
    ycycle=y(:,1:6);
    
    
    
%%%-----------------calculation of the PRC via integration of oscillators
%%%from unperturbed (y01) and perturbed (y02) initial points
    tphr=fi{i}*TT/2/pi;
    PPRC = [];
    for index =1:numvar
        for ii=1:length(fi{i})
            if tphr(ii)==0
                t=0;
                y=y0;
            else
                [t,y]=ode45(@mywillfun,[0 tphr(ii)],y0,odeset('Reltol',c1,'AbsTol',c2));
            end
            y01=y(end,:); %unperturbed initial point
            y02=y01+[zeros(1,index-1) dy zeros(1,length(x0)-index)]; %perturbed initial point
            tphr1=0;
            tphr2=0;
            
            for jj=1:3 %integration through 2-3 cycles
                [t1,y1]=ode45(@mywillfun,[0 10*tinitial],y01,odeset('events',@newcycle,'Reltol',c1,'AbsTol',c2));
                [t2,y2]=ode45(@mywillfun,[0 10*tinitial],y02,odeset('events',@newcycle,'Reltol',c1,'AbsTol',c2));
                y01=y1(end,:);
                y02=y2(end,:);
                tphr1=tphr1+t1(end);
                tphr2=tphr2+t2(end);
            end
            PPRC(ii,index)=(tphr1-tphr2)*2*pi/TT/dy; %calculation of phase response from time differences
            percon=2*pi/dy;
            PPRC=mod(PPRC+percon/2,percon)-percon/2; %correction for the case (tphr1-tphr2) accidentally contains a multiple of T
        end
    end
    
    % % %
    
    fi0 = ficycle;
    y0 = ycycle;
    
    y0=y0(:,1:6);
    fi{i}=fi{i}';
    
    figure(i);
    %if i==1
        title(fprintf('************FIGURE %d*************',i));
    %end
    subplot(2,1,1);
    
    hold on;
    plot(fi0,y0(:,1:3));
    title({'Time Series for E, L, and C Cells with a Tonic Drive of ' ve});
    xlabel('Phase (\phi)');
    legend('E Cell','L Cell','C Cell');
    hold off;
    subplot(2,1,2);
    plot(fi{i},PPRC);
    title({'Phase Response Curve with a Tonic Drive of ' ve});
    xlabel('Phase (\phi)');
    legend('E Cell','L Cell','C Cell');
    hold off;
    

    dfi=fi{i}(2)-fi{i}(1);
    dt=dfi/(2*pi/TT);
    y=interp1(fi0,y0,fi{i},'spline');
    shift{i}=-pi:pi/20:pi;%phase-shifts at which the coupling function is determined
    
    
    %connectivity matrix of segmental oscillators. Rows 1,2,...,6 mean
    %E,L,C,C,L,E cells, repectively; +1=excitatory connection; -1=inhibitory
    %connection
    conn=[0 0 0 -1 0 0; 1 0 0 -1 0 0;1 -1 0 -1 0 0];
    conn=[conn;0 0 0 0 0 0;1 0 0 0 0 0 ; 0 0 0 -1 0 0];
    conn=[conn;flipud(fliplr(conn))];
    
    % numerical integration of coupling function
    pert=[];
    HH=[];
    for kkk=1:length(shift{i}) %phase shift between coupled oscillators
        for iii=1:3 %the effect of incoming connections of cell i
            for j=1:length(fi{i})-1 %integrating over a period
                w=abs(conn(iii,:));%outgoing connection strengths of i^th cell
                
                fr=interp1(fi{i}',y,mod(fi{i}(j)+shift{i}(kkk),2*pi),'spline'); %calculating firing rates
                fr=max(fr,0);
                
                va=sign(conn(iii,:))-y(j,iii); %this term takes inzto account that connections have reversion potential (+1 or -1)
                pert(j,iii)=sum(fr.*w.*va*dt); %perturbations from incoming connections
            end
        end
        
    
    %The coupling function is determined from the integral of pert*PRC over
    %a period summed for all elements of the oscillator.
    h=pert.*PPRC(1:end-1,:);
    HH(kkk)=2*sum(sum(h)); %'2*' because only perturbations on half of the segment were integrated
    
    end
    
    H{i} = HH
    PRC{i} = PPRC
    T(i) = TT
    % % %     y0 =
    % % %     fi =
    % % %     PRC =
    % % %     T =
    % % %
    % % %
    % % %
    % % %
    % % %     shift{i} = % ????
    % % %
    % % %     H{i}=  % The coupling function is H(shift)
    % % %
    % % %     fi{i} = % ????
    % % %
    % % %     PRC{i} =  % Phase response curve is PRC(fi)
    % % %
    % % %     T(i) = % ????
    
    
    
    
    
end
figure(10);
for i=1:length(shift)
    plot(-shift{i},H{i});
    title('***********FIGURE 6*************');
    hold on;
end
legend('0.005','0.0075','0.01','0.02','0.04','0.06','0.07');

%%%%%%%%%%%%%%%%%stable zeros of coupling function
disp('calculating roots of copupling functions');
for i=1:length(shift)
    for j=1:length(H{i})-1
        if H{i}(j+1)*H{i}(j)<=0
            if H{i}(j+1)>H{i}(j)
                stab(i)=(shift{i}(j)*H{i}(j+1)-shift{i}(j+1)*H{i}(j))/(H{i}(j+1)-H{i}(j))/2/pi;
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%plotting Fig 7
figure(11);
plot(1./T,-stab,'s-');
title('***********FIGURE 7*************');
disp('Ready!');
