
npt=struct;
npt.RabiA = sqrt(2)*8.84;
npt.delta   = 200;
npt.d0      = npt.delta*2;
npt.eta     = .7;
npt.doFit   = 1;
out=threeLevelEvolve(npt);

hF = figure(10);
clf
hF.Color='w';
hF.Position=[100 100 500 250];
p1=plot([out.t],[out.rho11],'-','linewidth',1);
hold on
p2=plot([out.t],[out.rho22],'-','linewidth',1);
p3=plot([out.t],[out.rho33],'-','linewidth',1);
xlim([0 8]);
xlabel('time (ms)');
ylabel('population');
ylim([0 1]);

plot([0 8],[1 1]*out.EffectiveRabiAmplitude,'k-');
plot([1 1]*(1/out.EffectiveRabiTheory),[0 1]*out.EffectiveRabiAmplitude,'k-');

legend([p1 p2 p3],{'$\rho_{11}$','$\rho_{22}$','$\rho_{33}$'},'interpreter','latex');
set(gca,'fontsize',12,'xgrid','on','ygrid','on','box','on','linewidth',1);
str=['$\Omega = ' num2str(round(npt.RabiA,2)) '~\mathrm{kHz}$' newline ...
    '$\Delta=' num2str(npt.delta) '~\mathrm{kHz}$' newline ...
    '$\eta=' num2str(npt.eta) '$'];
text(.01,.98,str,'units','normalized','fontsize',12,'interpreter','latex','verticalalignment','top',...
    'backgroundcolor',[1 1 1 .5]);
%%


npt=struct;
npt.RabiA = sqrt(2)*8.84;
npt.delta   = 20;
npt.d0      = npt.delta*2;
npt.eta     = .4;
npt.doFit   = 1;
out=threeLevelEvolve(npt);

hF = figure(11);
clf
hF.Color='w';
hF.Position=[600 100 500 250];
p1=plot([out.t],[out.rho11],'-','linewidth',1);
hold on
p2=plot([out.t],[out.rho22],'-','linewidth',1);
p3=plot([out.t],[out.rho33],'-','linewidth',1);
xlim([0 1]);
xlabel('time (ms)');
ylabel('population');
ylim([0 1]);

plot([0 8],[1 1]*out.EffectiveRabiAmplitude,'k-');
plot([1 1]*(1/out.EffectiveRabiTheory),[0 1]*out.EffectiveRabiAmplitude,'k-');

legend([p1 p2 p3],{'$\rho_{11}$','$\rho_{22}$','$\rho_{33}$'},'interpreter','latex');
set(gca,'fontsize',12,'xgrid','on','ygrid','on','box','on','linewidth',1);
str=['$\Omega = ' num2str(round(npt.RabiA,2)) '~\mathrm{kHz}$' newline ...
    '$\Delta=' num2str(npt.delta) '~\mathrm{kHz}$' newline ...
    '$\eta=' num2str(npt.eta) '$'];
text(.01,.98,str,'units','normalized','fontsize',12,'interpreter','latex','verticalalignment','top',...
    'backgroundcolor',[1 1 1 .5]);
