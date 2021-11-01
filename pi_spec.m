function [t_rabi,rho_rabi] = pi_spec(npt)

if nargin~=1
   npt=struct;
   npt.mod_amp=2;
   npt.delta0=0;
   npt.Tp = 1/npt.mod_amp;
end

Tpi = 0.5*npt.Tp;

pi_rabi=@(t) npt.mod_amp*(t>-Tpi/2).*(t<Tpi/2);
pi_detune=npt.delta0;

rho0=[1 0 0 0];%start in excited state
%% Define equations of motion

function output=rabi_diffeq(t,y)
    M=(2*pi)*0.5*1j*[...
               0        0           pi_rabi(t)     -pi_rabi(t); ...
               0        0           -pi_rabi(t)    pi_rabi(t); ...
               pi_rabi(t)  -pi_rabi(t)    2*pi_detune 0  ; ...
              -pi_rabi(t)  pi_rabi(t)     0           -2*pi_detune];
    output=M*y;
end

%% Tevolve
tspan=npt.Tp*.75*[-1 1];

% Time evolve the chirp pulse
[t_rabi,rho_rabi]=ode45(@rabi_diffeq, tspan, rho0);

%% Plot Init
if npt.doPlot
    hf=figure(101);
    clf
    hf.Position = [800 200 400 400];
    hf.Color='w';

    axes
    co=get(gca,'colororder');
    hold on
    p1=plot(t_rabi,rho_rabi(:,1),'-','linewidth',1,'color',co(4,:),'linewidth',2);
    % plot([min(tspan) max(tspan)],[1 1]*P0,'k--');
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



