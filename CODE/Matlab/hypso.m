DEM = GRIDobj('G:/My Drive/INVESTIGACION/POSDOC/Data/Raster/subcuencas/dem_V5_CUENCA219_219.tif');
%DEM = inpaintnans(DEM);
FD = FLOWobj(DEM);

Hmax  = upslopestats(FD,DEM,'max');
Hmin  = upslopestats(FD,DEM,'min');
Hmean = upslopestats(FD,DEM,'mean');

S = STREAMobj(FD,'minarea',100);
z = getnal(S,DEM);
Hmean = getnal(S,Hmean);
Hmin  = getnal(S,Hmin);
Hmax  = getnal(S,Hmax);
HI = (Hmean - Hmin)./(Hmax - Hmin);

plotc(S,HI)
axis image
padextent(2000);
[x,y] = getoutline(DEM,true);
hold on
plot(x,y,'k')

St = trunk(klargestconncomps(S));
zt = getnal(St,DEM);
Hmeant = nal2nal(St,S,Hmean);
Hmaxt = nal2nal(St,S,Hmax);
Hmint = nal2nal(St,S,Hmin);
HIt   = nal2nal(St,S,HI);

figure
%plotdzshaded(St,[Hmax, Hmean],'facecolor',[0.3010 0.7450 0.9330])
hold on
%plotdzshaded(St,[Hmean Hmin],'facecolor',[0.9290 0.6940 0.1250]);
plotdz(St,DEM,'color','k');
yyaxis right
plotdz(St,HIt,'color','b','LineWidth',1.5)
ylabel('Hypsometric integral');
xlim([0 max(St.distance)])
legend('River profile', 'Hypsometric integral')
box on
