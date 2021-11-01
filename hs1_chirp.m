function [t_chirp,rho_chirp,t_hs1,rho_hs1]=hs1_chirp(nput)

if nargin~=1
   nput=struct;
   nput.Tp=5;     % Plse time in ms   
   nput.mod_amp=2;
   nput.freq_amp=5;
   nput.delta0=0;
end

% Beta Parameter
beta=asech(0.005);

% Amplitude of HS1
B=pi/(2*beta);  % Equal area
B=1;            % Equal amplitude

hs1_rabi=@(t) nput.mod_amp*sech(2*beta/nput.Tp*t);
hs1_detune=@(t) nput.freq_amp*tanh(2*beta/nput.Tp*t)+nput.delta0;

chirp_rabi=@(t) B*nput.mod_amp*(t>-nput.Tp/2).*(t<nput.Tp/2);
chirp_detune=@(t) -nput.freq_amp+2*((t+nput.Tp/2)/nput.Tp)*nput.freq_amp.*(t>-nput.Tp/2).*(t<nput.Tp/2)+2*nput.freq_amp.*(t>=nput.Tp/2)+nput.delta0;

rho0=[1 0 0 0];%start in excited state
%% Define equations of motion

function output=HS1_diffeq(t,y)
    M=(2*pi)*0.5*1j*[ 0  0  hs1_rabi(t)      -hs1_rabi(t); ...
               0  0 -hs1_rabi(t)    hs1_rabi(t); ...
               hs1_rabi(t) -hs1_rabi(t)  2*hs1_detune(t) 0  ; ...
              -hs1_rabi(t)  hs1_rabi(t)  0      -2*hs1_detune(t)];
    output=M*y;
end

function output=chip_diffeq(t,y)
    M=(2*pi)*0.5*1j*[ 0  0  chirp_rabi(t)      -chirp_rabi(t); ...
               0  0 -chirp_rabi(t)    chirp_rabi(t); ...
               chirp_rabi(t) -chirp_rabi(t)  2*chirp_detune(t) 0  ; ...
              -chirp_rabi(t)  chirp_rabi(t)  0      -2*chirp_detune(t)];
    output=M*y;
end



%% Tevolve

tspan=nput.Tp*.75*[-1 1];

opts_1 = odeset('RelTol',1e-2,'AbsTol',1e-3);

% Time evolve the chirp pulse
[t_chirp,rho_chirp]=ode45(@chip_diffeq, tspan, rho0,opts_1);


% Time evolve the HS1 pulse
[t_hs1,rho_hs1]=ode45(@HS1_diffeq, tspan, rho0);


% Landau-Zener prediction for chirp
f_rabi_chirp=B*nput.mod_amp;
ramp_ramp_chirp=2*nput.freq_amp/nput.Tp;
P0=1-exp(-pi^2*f_rabi_chirp^2/ramp_ramp_chirp);
%% Plot Init
if nput.doPlot

    hf=figure(100);
    clf
    hf.Position = [100 200 600 400];
    hf.Color='w';
    co=get(gca,'colororder');

    % Rabi frequency plot
    subplot(221)
    plot(t_chirp,hs1_rabi(t_chirp),'-','linewidth',2);
    hold on
    plot(t_chirp,chirp_rabi(t_chirp),'-','linewidth',2);
    set(gca,'fontsize',10,'box','on','linewidth',1);
    xlabel('time (ms)');
    ylabel('rabi frequency (kHz)');
    ylim([0 nput.mod_amp]);
    xlim([min(t_chirp) max(t_chirp)]);
    legend({'HS1','chirp'},'fontsize',8);



    % Detuning plot
    subplot(222)
    plot(t_chirp,hs1_detune(t_chirp),'-','linewidth',2);
    hold on
    plot(t_chirp,chirp_detune(t_chirp),'-','linewidth',2);
    set(gca,'fontsize',10,'box','on','linewidth',1);
    xlabel('time (ms)');
    ylabel('detuning (kHz)');
    ylim([-nput.freq_amp nput.freq_amp]*1.2+nput.delta0)
    xlim([min(t_chirp) max(t_chirp)]);
    legend({'HS1','chirp'},'fontsize',8,'location','southeast');


    subplot(223);
    hold on
    p1=plot(t_chirp,rho_chirp(:,1),'-','linewidth',1,'color',co(3,:),'linewidth',2);
    p2=plot(t_chirp,rho_chirp(:,2),'-','linewidth',1,'color',co(4,:),'linewidth',2);
    plot([min(tspan) max(tspan)],[1 1]*P0,'k--');
    ylabel('population');
    xlabel('time (ms)');
    xlim(tspan);
    ylim([0 1]);
    set(gca,'fontsize',10,'box','on');


    subplot(224);
    hold on
    p1=plot(t_hs1,rho_hs1(:,1),'-','linewidth',1,'color',co(3,:),'linewidth',2);
    p2=plot(t_hs1,rho_hs1(:,2),'-','linewidth',1,'color',co(4,:),'linewidth',2);
    ylabel('population');
    xlabel('time (ms)');
    xlim(tspan);
    ylim([0 1]);
    set(gca,'fontsize',10,'box','on');

end
%%
doSave=0;
if doSave
    fprintf('saving figures ...');
    
    % Save to png
    print(hf,'hs1_chirp_evolution.png','-dpng','-r400'); 

    disp('done');
end
end
