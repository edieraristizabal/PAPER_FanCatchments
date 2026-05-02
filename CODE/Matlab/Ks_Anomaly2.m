%% Load DEM and get stream network
DEM = GRIDobj('Nechi/Subcuencas/Aburra/porce_VdeA.tif');
FD  = FLOWobj(DEM);
A   = flowacc(FD);

%%
S   = STREAMobj(FD,A>5000);
S   = klargestconncomps(S);
%S   = trunk(S);

% Here we just take the standard value of theta (concavity)
theta = 0.54;

% Smoothing the river
K1    = 10; % for profile smoothing
K2    = 20; % for ksn smoothing
zs    = crs(S,DEM,'split',false,'K',K1);

% Uncertainties in river elevations 
s      = std(getnal(S,DEM) - zs);

% Ksn and mean Ksn
k  = ksn(S,zs,A,theta);
ks = smooth(S,k,'K',K2);

m  = mean(ks);
m  = getnal(S)+m;

% Calculate the model that Ksn has the same value everywhere
ghat  =  m.*(getnal(S,A)*DEM.cellsize^2).^(-theta);
zhat  =  cumtrapz(S,ghat);

%% plot
figure
plotdz(S,DEM)
hold on
plotdz(S,zhat)
hold off

%% Simulation. Note that I use parfor here because parallelization
% speeds up the computation tremendously.
nsim   = 100;
ksnsim = cell(1,nsim);
parfor r = 1:nsim
    zsim = zhat + randn(numel(zs),1)*s; 
    zsim = crs(S,zsim,'split',false,'K',K1);
    ksim = ksn(S,zsim,A,theta);
    ksnsim{r} = smooth(S,ksim,'K',K2);
end
hold on
ksnsim  = horzcat(ksnsim{:});
ksnperc = prctile(ksnsim,[5 95],2);

%% Finally, plot it
figure
plotdz(S,DEM);
yyaxis right
plotdzshaded(S,[max(ks,m) m],'FaceColor','r');
hold on
plotdzshaded(S,[min(ks,m) m],'FaceColor','b');
plotdz(S,ks,'color',[.6 .6 .6])
plotdzshaded(S,ksnperc,'facecolor',[.5 .5 .5]);
plotdz(S,m,'color','k')
ylabel('K_{sn}')
xlim([0 max(S.distance)])