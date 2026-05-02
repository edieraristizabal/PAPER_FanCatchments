%% Load DEM and get stream network
DEM = GRIDobj('../../Nechi/TIF/porce_VdeA.tif');
FD  = FLOWobj(DEM);
A   = flowacc(FD);

%%
S   = STREAMobj(FD,A>10000);
S   = klargestconncomps(S);

%%
Sm = modify(S,'streamorder','<5');

%%
Sm = modify(S,'distance',[50000 100000]);

%%
Sm = modify(S,'upstreamto',DEM>1500);

%%
Sm = modify(S,'tributaryto',trunk(S));

%%
Sm = modify(S,'interactive','reachselect');
plot(Sm)