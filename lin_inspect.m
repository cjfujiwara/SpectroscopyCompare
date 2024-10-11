% Simulation Settings

npt=struct;             % Initialize the structure
npt.mod_amp     = 8;    % [kHz] Modulation Rabi Amplitude
npt.freq_amp    = 25;    % [kHz] Frequency Detuning Amplitude
npt.delta0      = 0;     % [kHz] Center Frequency
npt.doPlot      = 0;    % Show the time traces?
npt.Tp          = 16.6;     % [ms] Pulse time
npt.LinRampTime = .5;

% How many detunings to simulate
tVec = linspace(0,1,100);
tVec(1)=[];

% Initialize density vectors for chirp, hs1, pi
y_chirp=zeros(length(tVec),2);

%% Evlove TDSE for each detuning

for kk=1:length(tVec)
    fprintf(['(' num2str(kk) ' of ' num2str(length(tVec)) ') ' ...
        'T = ' num2str(tVec(kk)) ' ... ']);
    npt.Tp=tVec(kk);
    [t_chirp,rho_chirp,t_hs1,rho_hs1]=hs1_chirp(npt);
    y_chirp(kk,:)=rho_chirp(end,1:2);
   disp('done');
end

%% Plot the Results
str=['$\Omega_0=2\pi \times' num2str(npt.mod_amp) '~\mathrm{kHz}$' newline ...
    '$\Delta_{\mathrm{amp}}=' num2str(npt.freq_amp) '~\mathrm{kHz}$'];


hf1=figure;
clf
hf1.Color='w';
hf1.Position(3:4)=[400 400];

co=get(gca,'colororder');
plot(tVec,y_chirp(:,2),'linewidth',1,'color',co(1,:),'linewidth',1);

xlabel('sweep time (ms)');
ylabel('population transfer');
legend({'linear chirp'});

text(.01,.98,str,'units','normalized','interpreter','latex',...
    'verticalalignment','top');
ylim([0 1.2]);
