function out = threeLevel(delta)

npt=struct;
npt.RabiA = sqrt(2)*8.84;

if nargin~=1
    delta = 33;
end

etaVec = 3:.1:.3;

etaVec=.3;

for kk=1:length(etaVec)
    npt.delta   = delta;
    npt.d0      = npt.delta*2;
    npt.eta     = etaVec(kk);
    npt.doFit   = 1;
    out(kk)=threeLevelEvolve(npt);
end

end

