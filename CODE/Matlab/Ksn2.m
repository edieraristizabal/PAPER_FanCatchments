%DEM = GRIDobj('../../Data/AOI_60m_full.tif');
DEM = fillsinks(DEM)
FD  = FLOWobj(DEM,'preprocess','carve');
A   = flowacc(FD);

%%
S   = STREAMobj(FD,A > 1000);
S = modify(S,'distance',100);

%% channel network
S = klargestconncomps(S,1);
imageschs(DEM);
hold on
plot(S,'k')
hold off

%% longitudinal profile
plotdz(S,DEM)

%% calculate
STATS = slopearea(S,DEM,A)

G   = gradient8(DEM);
KSN = G./(A.*(A.cellsize^2)).^STATS.theta;

[x,y,ksn] = STREAMobj2XY(S,KSN);
scatter(x,y,5,ksn,'filled')
caxis([0 50]) % ajustar escala
title('ksn values')
box on
colorbar
axis image

%%

MS = STREAMobj2mapstruct(S,'seglength',1000,'attributes',...
   {'ksn' KSN @mean ...
    'uparea' (A.*(A.cellsize^2)) @mean ...
    'gradient' G @mean});

symbolspec = makesymbolspec('line',...
    {'ksn' [min([MS.ksn]) max([MS.ksn])] 'color' jet(6)});

imageschs(DEM,DEM,'colormap',gray,'colorbar',false);
mapshow(MS,'SymbolSpec',symbolspec);
colorbar
axis image

shapewrite(MS,'prueba.shp');