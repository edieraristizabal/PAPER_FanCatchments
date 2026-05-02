%% map viwer
%DEM = GRIDobj('../../Nechi/nechi_30m.tif');
%FD = FLOWobj(DEM);
S  = STREAMobj(FD,'minarea',20000);
mappingapp(DEM,S)