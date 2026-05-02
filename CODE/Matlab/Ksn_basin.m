DEM = GRIDobj('../../Nechi/TIF/porce_VdeA.tif');
DEM = inpaintnans(DEM);
FD  = FLOWobj(DEM);
DEM = imposemin(FD,DEM);
A   = flowacc(FD);
S   = STREAMobj(FD,A > 1000);

%% calculate Ksn
k   = ksn(S,DEM,A,0.67);

%% to assign basins
D   = drainagebasins(FD,S);
d   = getnal(S,D);
km  = accumarray(d,k,[],@mean);

%% plot
K   = GRIDobj(DEM)*nan;
K.Z(D.Z>0) = km(D.Z(D.Z>0));
imageschs(DEM,K,'colorbarylabel','K_{sn} [m^0.9]')