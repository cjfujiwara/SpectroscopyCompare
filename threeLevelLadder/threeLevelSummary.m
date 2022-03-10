deltaVec = [10 20 30 40 50 100 200];



for kk=1:length(deltaVec)
    out{kk} = threeLevel(deltaVec(kk));
end

%%
hF = figure(3);
clf
hF.Color='w';
hF.Position=[100 100 600 400];
for kk=1:length(out)
    plot([out{kk}.eta],[out{kk}.EffectiveRabiFit]./[out{kk}.EffectiveRabiTheory],'linewidth',2);
    hold on
end

xlabel('\eta');
ylabel('$\tilde{\Omega}_{\mathrm{eff}}$ numeric/$\tilde{\Omega}_{\mathrm{eff}}$ theory','interpreter','latex')

legStr={'$\Delta = 10~\mathrm{kHz}$','$\Delta = 20~\mathrm{kHz}$','$\Delta = 30~\mathrm{kHz}$','$\Delta = 40~\mathrm{kHz}$','$\Delta = 50~\mathrm{kHz}$','$\Delta = 100~\mathrm{kHz}$','$\Delta = 200~\mathrm{kHz}$'};

legend(legStr,'location','eastoutside','interpreter','latex','fontsize',12);
set(gca,'fontsize',12,'xgrid','on','ygrid','on');
xlim([.1 1]);