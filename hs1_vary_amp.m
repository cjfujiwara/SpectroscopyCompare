% Simulation Settings

npt=struct;             % Initialize the structure
npt.Tp          = 16.6;    % [ms] Pulse time
npt.freq_amp    = 25;  % [kHz] Frequency Detuning Amplitude
npt.delta0      = 0;    % [kHz] Center Frequency
npt.doPlot      = 0;    % Show the time traces?

% How many detunings to simulate
modamp_vec=linspace(.1,10,100);

% Initialize density vectors for chirp, hs1, pi
y_chirp=zeros(length(delta0vec),2);
y_hs1=zeros(length(delta0vec),2);
y_pi=zeros(length(delta0vec),2);

%% Evlove TDSE for each detuning

for kk=1:length(modamp_vec)
    fprintf(['(' num2str(kk) ' of ' num2str(length(modamp_vec)) ') ' ...
        'modamp = ' num2str(modamp_vec(kk)) ' ... ']);
    npt.delta0=0;
    npt.mod_amp = modamp_vec(kk);
    
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
str=['$\Delta_0=2\pi \times' num2str(npt.delta0) '~\mathrm{kHz}$' newline ...
    '$T_p=' num2str(npt.Tp) '~\mathrm{ms}$' newline ...
    '$\Delta_{\mathrm{amp}}=' num2str(npt.freq_amp) '~\mathrm{kHz}$'];


hf1=figure;
clf
hf1.Color='w';
hf1.Position(3:4)=[400 400];

co=get(gca,'colororder');
% plot(modamp_vec,y_chirp(:,2),'linewidth',1,'color',co(1,:),'linewidth',1);
hold on
plot(modamp_vec,y_hs1(:,2),'linewidth',1,'color',co(2,:),'linewidth',1);
% plot(modamp_vec,y_pi(:,2),'linewidth',1,'color',co(4,:),'linewidth',1);

xlabel('peak rabi (kHz)');
ylabel('population transfer');
% legend({'linear chirp','HS1','\pi'});
legend('HS1','location','southeast');

text(.01,.98,str,'units','normalized','interpreter','latex',...
    'verticalalignment','top');

set(gca,'box','on','linewidth',1)

