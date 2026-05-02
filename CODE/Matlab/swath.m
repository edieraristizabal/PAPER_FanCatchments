%DEM = GRIDobj('../../Nechi/TIF/porce_VdeA.tif');
%FD = FLOWobj(DEM);
S  = STREAMobj(FD,'minarea',1e6,'unit','map');
S = klargestconncomps(S);
S = trunk(S);
%%
imageschs(DEM,[],'colormap',[1 1 1],'colorbar',false);
hold on
plot(S)

%%
DEM = inpaintnans(DEM);
z = getnal(S,DEM);
[zmax,mask] = maplateral(S,DEM,2000,@max);

subplot(2,1,1)
imageschs(DEM,mask,'truecolor',[1 0 0],...
                   'colorbar',false,...
                   'ticklabels','nice');

subplot(2,1,2)               
plotdz(S,z)
hold on
plotdz(S,zmax)
legend('River profile','maximum heights in 2 km distance')

figure
zmaxs = crs(S,zmax,'mingradient',nan,'K',5,'tau',0.99,'split',0);
plotdz(S,z); 
hold on
plotdzshaded(S,[zmaxs z],'FaceColor',[.6 .6 .6])
hold off

%%
figure 
SW = STREAMobj2SWATHobj(S,DEM,'width',4000);
plotdz(SW)
xlabel('Distance upstream [m]');
ylabel('Elevation [m]');