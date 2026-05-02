%% Section 1 import dem
DEM = GRIDobj('G:/My Drive/INVESTIGACION/POSDOC/Data/Raster/subcuencas/dem_V5_CUENCA292_292.tif')
%DEM.Z(DEM.Z<-9998)=NaN; 

% Section 2 Flow direction, stream network, drainage area
FD = FLOWobj(DEM,'preprocess','carve');
A = flowacc(FD);
S = STREAMobj(FD,A>100);
S   = removeedgeeffects(S,FD);
D  = drainagebasins(FD,S);

%cn  = griddedcontour(DEM,[500 500],true);
%S  = modify(S,'upstreamto',cn);
%s = trunk(klargestconncomps(S)); %main drainage

% plot the channel network
imageschs(DEM,[],'colormap',[1 1 1],'colorbar',false);
hold on
plot(S)

% calculate chi
c=chitransform(S,A,'mn',0.24,'a0',1e6);

% plot chi map
imageschs(DEM,[],'colormap',[1 1 1],'colorbar',false);
hold on
plotc(S,c)
%%
% export table
%t = table;
%[t.lat,t.lon,t.c] = STREAMobj2latlon(S,c);
%writetable(t,'Nechi/Subcuencas/Nechi/Chi_nechi_mn024.xlsx');


