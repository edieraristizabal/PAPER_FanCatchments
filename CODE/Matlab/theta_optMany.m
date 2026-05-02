%code to find the best concavity value for MANY BASINS
clear
clc
DEM = GRIDobj('G:/My Drive/INVESTIGACION/POSDOC/Data/Raster/subcuencas/dem_V5_CUENCA292_292.tif');
FD  = FLOWobj(DEM);
A   = flowacc(FD);
S   = STREAMobj(FD,A>1000);
%S   = klargestconncomps(S);
plotdz(S,DEM)
%
mn = mnoptim(S,DEM,A,'mnrange',[0.1 1],'optvar','mn','crossval',false,'plot',false);
mn1 = mn.mn
mn2 = mnoptimvar(S,DEM,A,'varfun',@robustcov)

subplot(2,1,1)
c1 = chitransform(S,A,'mn',mn1);
plotdz(S,DEM,'distance',c1);
h(1) = subplotlabel(gca,['mn = ' num2str(mn1)]);
subplot(2,1,2)
c2 = chitransform(S,A,'mn',mn2);
plotdz(S,DEM,'distance',c2);
h(2) = subplotlabel(gca,['mn = ' num2str(mn2)]);
bigger(h,6)
