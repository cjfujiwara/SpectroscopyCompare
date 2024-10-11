function out = threeLevelEvolve(npt)

% Eta Parameter
eta=npt.eta;

% Define time depdendent pulses
d   = @(t) npt.delta;
A   = @(t) npt.RabiA;
B   = @(t) npt.eta*npt.RabiA;
d0  = npt.d0;

Cbare=A(0)*B(0)/(2*d(0));

deff = (2*d(0)-d0)+(A(0)^2-B(0)^2)/(4*d(0));

Ceff = sqrt(Cbare.^2+deff.^2);
Tc = 1/Ceff;

disp(['effective rabi '  num2str(Cbare)]);
disp(['effective detuning ' num2str(deff)]);
disp(['bare period ' num2str(1/Cbare)]);
disp(['detuning effective rabi ' num2str(Ceff)]);
disp(['detuning effective period ' num2str(Tc)]);

% Define time
t1 = 0;
if ~isnan(Tc)
    t2 = 4*Tc;
else
    t2 = 3/A(0);
end
tspan = [t1 t2];

% Intial Density matrix
rho0 = [1 0 0;
    0 0 0;
    0 0 0];

function output=diffeq(t,rhoVec)    
    % The Hamiltonian at this time
    H = 2*pi*0.5*[...
        0       A(t)    0;
        A(t)    -2*d(t) B(t);
        0       B(t)    -4*d(t)+2*d0];    
    
    % Convert density vector in density matrix
    rhoMat = reshape(rhoVec,[3 3]);
    
    % Time Evolve according to liouville von nueman 
    drhoMatdt = -1i*(H*rhoMat-rhoMat*H);
    
    % Convert back into vector
    output=drhoMatdt(:);
end


opts_1 = odeset('RelTol',1e-4,'AbsTol',1e-5,'MaxStep',.5e-3);

% Time evolve the chirp pulse
[t_out,rho_out]=ode45(@diffeq, tspan, rho0(:),opts_1);

hF=figure(1);
hF.Color='w';
clf
co=get(gca,'colororder');

subplot(3,3,[1 2]);
plot(t_out,real(rho_out(:,1)),'-','color',co(1,:),'linewidth',2);
xlabel('time (ms)');
ylabel('\rho_{00}');
xlim(tspan);
set(gca,'box','on','linewidth',1,'fontsize',10,'xgrid','on','ygrid','on');
ylim([0 1]);

subplot(3,3,[4 5]);
plot(t_out,real(rho_out(:,5)),'-','color',co(2,:),'linewidth',2);
xlabel('time (ms)');
ylabel('\rho_{11}');
xlim(tspan);
set(gca,'box','on','linewidth',1,'fontsize',10,'xgrid','on','ygrid','on');
ylim([0 1]);

ax=subplot(3,3,[7 8]);
plot(t_out,real(rho_out(:,9)),'-','color',co(3,:),'linewidth',2);


xlabel('time (ms)');
ylabel('\rho_{22}');
xlim(tspan);
set(gca,'box','on','linewidth',1,'fontsize',10,'xgrid','on','ygrid','on');
ylim([0 1.2]*sqrt(Cbare^2/(deff^2+Ceff^2)));

hold on
plot([1 1]*Tc,[0 1],'k--');
plot([0 1]*3*Tc,[1 1]*sqrt(Cbare^2/(deff^2+Ceff^2)),'k--');


ax3=subplot(1,3,3);
ax3.Units='pixels';
pp=ax3.Position;
delete(ax3);

t1=uitable('position',pp,'units','pixels','rowname',{},'columnname',{},...
    'ColumnWidth',{140 50});

t1.Data{1,1}='single detuning (kHz)';
t1.Data{2,1}='Bare Rabi (kHz)';
t1.Data{3,1}='eta';
t1.Data{4,1}='third state detuning (kHz)';

t1.Data{1,2}=num2str(d(0));
t1.Data{2,2}=num2str(A(0));
t1.Data{3,2}=num2str(eta);
t1.Data{4,2}=num2str(d0);

t1.Data{6,1}='bare two rabi (kHz)';
t1.Data{7,1}='effective detuning (kHz)';
t1.Data{8,1}='effective rabi (kHz)';
t1.Data{9,1}='period (ms)';
t1.Data{10,1}='amplitude (arb)';

t1.Data{6,2}=num2str(Cbare);
t1.Data{7,2}=num2str(deff);
t1.Data{8,2}=num2str(Ceff);
t1.Data{9,2}=num2str(Tc);
t1.Data{10,2}=num2str(sqrt(Cbare^2/(deff^2+Ceff^2)));

out = struct;
out.eta = eta;
out.RabiBare = npt.RabiA;
out.Delta = npt.delta;
out.Delta0 = npt.d0;

out.t = t_out;
out.rho11 = real(rho_out(:,1));
out.rho22 = real(rho_out(:,5));
out.rho33 = real(rho_out(:,9));

out.EffectiveRabiTheory = Ceff;
out.EffectiveRabiAmplitude = sqrt(Cbare^2/(deff^2+Ceff^2));
out.EffectiveRabiTheory2 = (sqrt((out.Delta)^2+(1+npt.eta^2)*npt.RabiA^2)-(out.Delta))/2;

plot([1 1]/out.EffectiveRabiTheory2,[0 1],'b--');

if npt.doFit
   myfit = fittype('A*sin(pi*f*t).^2','independent','t','coefficients',{'A','f'});
   fitopt=fitoptions(myfit);
   
   D2 = (2*out.Delta-out.Delta0);
    disp(D2-out.Delta)
f=out.EffectiveRabiTheory2;
f=sqrt(f^2+D2^2)+D2*((1-eta)+0.5*(1-eta)^2);
   
    plot([1 1]/f,[0 1],'g--');
   t1.Data{14,2}=num2str(f);


   fitopt.StartPoint=[max(real(rho_out(:,9))) f];
   fitopt.Upper=fitopt.StartPoint*1.5;
   fitopt.Low=fitopt.StartPoint.*[.9 .5];

   t1.Data{13,1}='freq theory2 (kHz)';
   t1.Data{13,2}=num2str(out.EffectiveRabiTheory2);
  
   fout = fit(t_out,real(rho_out(:,9)),myfit,fitopt);
   axes(ax)
   hold on
   plot(fout)
   CeffFit = fout.f;
   
   t1.Data{12,1}='freq fit (kHz)';
   t1.Data{12,2}=num2str(CeffFit);
   drawnow;
   out.EffectiveRabiFit = CeffFit;

end
% 
% figure(2)
% clf
% plot(t_out,(out.rho33-out.rho11)-out.rho22)

end

