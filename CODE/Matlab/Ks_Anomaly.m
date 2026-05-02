%% Load DEM
DEM = GRIDobj('../../Nechi/TIF/porce_VdeA.tif');
FD  = FLOWobj(DEM);
A = flowacc(FD);

%%
S   = STREAMobj(FD,A > 1000);
S = klargestconncomps(S);
%S = trunk(S);

%% Smooth elevations along the river using crs
%zs = crs(S,DEM,'split',false,'K',10);
% Calculate ksn
k = ksn(S,DEM,A);
% And smooth the output again a little. 
%ks = smooth(S,k,'K',20);

%% Plot the results
% First, plot the profile
plotdz(S,DEM);
%Add a axis on the right side
yyaxis right
% Calculate mean ksn
m = mean(k);
% Calculate a node-attribute list of mean Ksn
m = getnal(S)+m;
% Plot positive deviations from the mean
plotdzshaded(S,[max(k,m) m],'FaceColor','r');
hold on
% ... and negative ones
plotdzshaded(S,[min(k,m) m],'FaceColor','b');
% Plot some outlines because it looks nicer
plotdz(S,k,'color',[.6 .6 .6])
% Emphasize the mean value
plotdz(S,m,'color','k')
hold off
ylabel('K_{sn}')
xlim([0 max(S.distance)])