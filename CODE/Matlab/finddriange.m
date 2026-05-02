
% import dem and create flow direction and drainage network and export data
DEM = GRIDobj('G:/My Drive/INVESTIGACION/POSDOC/Data/Raster/subcuencas/dem_V5_CUENCA219_219.tif');
DEM = inpaintnans(DEM);
FD = FLOWobj(DEM);
A   = flowacc(FD);
S = STREAMobj(FD,A>0);
zt = repmat(getnal(S),1,3);

L   = readgeotable('G:\My Drive\INVESTIGACION\POSDOC\Data\Vector\recents.geojson');

L.Shape

zt = maplateral(L,DEM)