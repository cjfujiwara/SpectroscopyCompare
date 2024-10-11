% Simulation Settings

npt=struct;             % Initialize the structure
npt.mod_amp     = 8;    % [kHz] Modulation Rabi Amplitude
npt.freq_amp    = 20;     % [kHz] Frequency Detuning Amplitude
npt.delta0      = 0;     % [kHz] Center Frequency
npt.doPlot      = 0;    % Show the time traces?
% npt.Tp          = 2*npt.freq_amp/10;     % [ms] Pulse time
npt.Tp          = 1;     % [ms] Pulse time

npt.LinRampTime = .1;
 % npt.Tp = npt.freq_amp/10;

% How many detunings to simulate
delta0vec=linspace(-4*npt.freq_amp,4*npt.freq_amp,200);

% Initialize density vectors for chirp, hs1, pi
y_chirp=zeros(length(delta0vec),2);
y_hs1=zeros(length(delta0vec),2);
y_pi=zeros(length(delta0vec),2);

%% Evlove TDSE for each detuning

for kk=1:length(delta0vec)
    fprintf(['(' num2str(kk) ' of ' num2str(length(delta0vec)) ') ' ...
        'delta = ' num2str(delta0vec(kk)) ' ... ']);
    npt.delta0=delta0vec(kk);
    [t_chirp,rho_chirp,t_hs1,rho_hs1]=hs1_chirp(npt);
    y_chirp(kk,:)=rho_chirp(end,1:2);
    y_hs1(kk,:)=rho_hs1(end,1:2);   
    
    rabi_npt = npt;
    rabi_npt.Tp = 1/rabi_npt.mod_amp;
    [t_pi,rho_pi] = pi_spec(rabi_npt);
    y_pi(kk,:)=rho_pi(end,1:2);   
    disp('done');
end

%% Plot the Results
str=['$\Omega_0=2\pi \times' num2str(npt.mod_amp) '~\mathrm{kHz}$' newline ...
    '$T_p=' num2str(npt.Tp) '~\mathrm{ms}$' newline ...
    '$\Delta_{\mathrm{amp}}=' num2str(npt.freq_amp) '~\mathrm{kHz}$' newline ...
    '$T_{\mathrm{lin ramp}} = ' num2str(npt.LinRampTime) '~\mathrm{ms}$'];


hf1=figure;
clf
hf1.Color='w';
hf1.Position(3:4)=[400 400];

co=get(gca,'colororder');
plot(delta0vec,y_chirp(:,2),'linewidth',1,'color',co(1,:),'linewidth',1);
hold on
plot(delta0vec,y_hs1(:,2),'linewidth',1,'color',co(2,:),'linewidth',1);
plot(delta0vec,y_pi(:,2),'linewidth',1,'color',co(4,:),'linewidth',1);

xlabel('detuning (kHz)');
ylabel('population transfer');
legend({'linear chirp','HS1','\pi'});

text(.01,.98,str,'units','normalized','interpreter','latex',...
    'verticalalignment','top');

