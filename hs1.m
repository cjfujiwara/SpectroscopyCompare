function hs1

Tp=10;
beta=asech(.005);
mod_amp=2;
freq_amp=5;
tt=linspace(-Tp*.75,Tp*.75,1e3);

B=pi/(2*beta);
B=1;
hs1_rabi=@(t) mod_amp*sech(2*beta/Tp*t);
hs1_detune=@(t) freq_amp*tanh(2*beta/Tp*t);

chirp_rabi=@(t) B*mod_amp*(t>-Tp/2).*(t<Tp/2);
chirp_detune=@(t) -freq_amp+2*((t+Tp/2)/Tp)*freq_amp.*(t>-Tp/2).*(t<Tp/2)+2*freq_amp.*(t>=Tp/2);

str_hs1_rabi='$\Omega_0~\mathrm{sech}\left[2\beta\frac{t}{T_p}\right]$';
str_chirp_rabi='$\Omega_0$';

str_hs1_detune='$f_0~\mathrm{tanh}\left[2\beta\frac{t}{T_p}\right]$';
str_chirp_detune='$f_0$';


%% Plot of power and detunings

hF1=figure(1);
clf
hF1.Color='w';
hF1.Position(3:4)=[400 400];

subplot(211)
plot(tt,hs1_rabi(tt),'-','linewidth',2);
hold on
plot(tt,chirp_rabi(tt),'-','linewidth',2);
set(gca,'fontsize',10,'box','on','linewidth',1);
xlabel('time');
ylabel('rabi frequency');
ylim([0 mod_amp*1.4]);
xlim([min(tt) max(tt)]);
legend({str_hs1_rabi,str_chirp_rabi},'interpreter','latex','fontsize',10,...
    'orientation','horizontal');

subplot(212)
plot(tt,hs1_detune(tt),'-','linewidth',2);
hold on
plot(tt,chirp_detune(tt),'-','linewidth',2);
set(gca,'fontsize',10,'box','on','linewidth',1);
xlabel('time');
ylabel('detuning');
legend({str_hs1_detune,str_chirp_detune},'interpreter','latex','fontsize',10,...
    'orientation','vertical','location','southeast');
ylim([-freq_amp freq_amp]*1.2)
xlim([min(tt) max(tt)]);

%%
tt=linspace(-3*Tp,3*Tp,1e6);
dt=tt(2)-tt(1);
L=length(tt);
Fs=1/dt;

f = Fs*(0:(L/2))/L;

fc=1E3; % Carrier frequency

hs1_yy=hs1_rabi(tt).*sin(2*pi*(fc+hs1_detune(tt)).*tt);
chirp_yy=chirp_rabi(tt).*sin(2*pi*(fc+chirp_detune(tt)).*tt);

% Testing FFT
S = 0.7*sin(2*pi*50*tt) + sin(2*pi*120*tt);
X = S + 2*randn(size(tt));
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

HS1_Y=fft(hs1_yy);
HS1_P2=abs(HS1_Y/L);
HS1_P1=HS1_P2(1:L/2+1);
HS1_P1(2:end-1) = 2*HS1_P1(2:end-1);

CHIRP_Y=fft(chirp_yy);
CHIRP_P2=abs(CHIRP_Y/L);
CHIRP_P1=CHIRP_P2(1:L/2+1);
CHIRP_P1(2:end-1) = 2*CHIRP_P1(2:end-1);

hF2=figure(2);
hF2.Position(3:4)=[400 400];
clf
hF2.Color='w';
plot(f-fc,HS1_P1,'linewidth',1);
% plot(f,P1);

xlim([-1 1]*freq_amp*5);

hold on
plot(f-fc,CHIRP_P1,'linewidth',1);


xlabel('detuning (kHz)');
ylabel('amplitude (au)');

legend({'HS1','Chip'});
%%
doSave=0;
if doSave
    fprintf('saving figures ...');
    
    % Save to png
    print(hF1,'hs1_chirp_time.png','-dpng','-r400'); 
    print(hF2,'hs1_chirp_fft.png','-dpng','-r400'); 

    disp('done');
end
end

