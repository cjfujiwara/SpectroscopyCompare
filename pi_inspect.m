% Simulation Settings

npt=struct;             % Initialize the structure
npt.mod_amp     = 8;    % [kHz] Modulation Rabi Amplitude
npt.delta0      = 0;     % [kHz] Center Frequency
npt.doPlot      = 0;    % Show the time traces?

% How many detunings to simulate
delta0vec=linspace(-5*npt.mod_amp,5*npt.mod_amp,1000);

% Initialize density vectors for chirp, hs1, pi

y_pi=zeros(length(delta0vec),2);

%% Evlove TDSE for each detuning

for kk=1:length(delta0vec) 
        fprintf(['(' num2str(kk) ' of ' num2str(length(delta0vec)) ') ' ...
        'delta = ' num2str(delta0vec(kk)) ' ... ']);
    npt.delta0=delta0vec(kk);

    rabi_npt = npt;
    rabi_npt.Tp = 1/rabi_npt.mod_amp;
    rabi_npt.mod_amp=8;
    [t_pi,rho_pi] = pi_spec(rabi_npt);
    y_pi(kk,:)=rho_pi(end,1:2);   
    disp('done');
end

%% Plot the Results
str=['$\Omega_0=2\pi \times' num2str(rabi_npt.mod_amp) '~\mathrm{kHz}$' newline ...
    '$T_p=' num2str(rabi_npt.Tp/2) '~\mathrm{ms}$'];


hf1=figure;
clf
hf1.Color='w';
hf1.Position(3:4)=[400 400];

co=get(gca,'colororder');

plot(delta0vec,y_pi(:,2),'linewidth',1,'color',co(4,:),'linewidth',1);

xlabel('detuning (kHz)');
ylabel('population transfer');
legend({'\pi'});


text(.01,.98,str,'units','normalized','interpreter','latex',...
    'verticalalignment','top');
