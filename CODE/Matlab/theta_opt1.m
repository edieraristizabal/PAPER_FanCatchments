% Code to find concavity best value for ONLY ONE BASIN
%DEM = GRIDobj('../../Nechi/TIF/esperanza_12m.tif');
%inpaintnans(DEM)
%FD = FLOWobj(DEM);
A  = flowacc(FD);
S  = STREAMobj(FD,A>10000);
%C  = griddedcontour(DEM,[500 500],true);
%S  = modify(S,'upstreamto',C);

[mn,results] = mnoptim(S,DEM,A,'mnrange',[0.1 1],'optvar','mn','crossval',false);
results = resume(results);
bestPoint(results)