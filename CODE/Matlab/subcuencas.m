% code to get morphometric parameters by Edier Aristizabal (2023)

% import dem and create flow direction and drainage network and export data
DEM = GRIDobj('G:/My Drive/INVESTIGACION/POSDOC/Nare/Raster/rionegro.tif');
FD = FLOWobj(DEM);
A   = flowacc(FD);
S = STREAMobj(FD,A>0); %cuando no se trabaja con A>0 se pierden muchos movimientos en masa pq ya no les ubica un drenaje
Sp = trunk(S);

% to calculate S-A plot data and export
a = getnal(S,A)*(DEM.cellsize^2);
g = gradient(S,DEM);
t1 = table;
[t1.flowacc] = a;
[t1.gradient] = g;

%writetable(t1,'G:/My Drive/INVESTIGACION/POSDOC/TopoToolbox/subcuencas/cartama_12m_SA.csv','WriteMode','overwrite');

% to export table with longitudinal profile network 
z = getnal(S,DEM);
t2 = table;
[t2.distance] = S.distance;
[t2.elevation] = z;

%recent
Lrec = GRIDobj('G:/My Drive/INVESTIGACION/POSDOC/Nare/Raster/rionegro_recent2.tif');
Lrec_h = getnal(S,Lrec).*getnal(S,DEM);
Lrec_d = getnal(S,Lrec).*S.distance;
[t2.Lrec_d] = Lrec_d;
[t2.Lrec_h] = Lrec_h;

%relict
%Lrel = GRIDobj('G:/My Drive/INVESTIGACION/POSDOC/Nechi/Subcuencas/Aburra/Raster/iguana_relict.tif');
%Lrelh = getnal(S,Lrel).*getnal(S,DEM);
%Lreld = getnal(S,Lrel).*S.distance;
%[t2.Lreld] = Lreld;
%[t2.Lrelh] = Lrelh;

writetable(t2,'G:/My Drive/INVESTIGACION/POSDOC/TopoToolbox/subcuencas/rionegro_network.csv','WriteMode','overwrite');
%%
% estimate HI for the drinage network
Hmax  = upslopestats(FD,DEM,'max');
Hmin  = upslopestats(FD,DEM,'min');
Hmean = upslopestats(FD,DEM,'mean');
Hmean = getnal(S,Hmean);
Hmin  = getnal(S,Hmin);
Hmax  = getnal(S,Hmax);
HI = (Hmean - Hmin)./(Hmax - Hmin);

% to obtain main drinage and its HI values and export
Sm = trunk(klargestconncomps(S));
%Sm = modify(S,'interactive','reachselect');

zt = getnal(Sm,DEM);
gt = gradient(Sm,DEM);
Hmeant = nal2nal(Sm,S,Hmean);
Hmaxt = nal2nal(Sm,S,Hmax);
Hmint = nal2nal(Sm,S,Hmin);
HIt   = nal2nal(Sm,S,HI);
t3 = table;
[t3.distance] = Sm.distance;
[t3.elevation] = zt;
[t3.gradient] = gt;
[t3.hypso] = HIt;

writetable(t3,'G:/My Drive/INVESTIGACION/POSDOC/TopoToolbox/subcuencas/rionegro_12m_main.csv','WriteMode','overwrite');

% to calculate hypsometric curve and export
[rf,elev]=hypscurve(DEM);
t4 = table;
[t4.hyp_elev] = elev;
[t4.hyp_rf] = rf;

writetable(t4,'G:/My Drive/INVESTIGACION/POSDOC/TopoToolbox/subcuencas/rionegro_12m_hypso.csv','WriteMode','overwrite');

% Calculate knickpoints and export
[~,kp] = knickpointfinder(Sp,DEM,'tol',100,'split',false,'plot',true); % el valor representa la altura minima en m de lso knickpoints
t5 = table;
[t5.x] = kp.x;
[t5.y] = kp.y;
[t5.z] = kp.z;
[t5.dz]= kp.dz;
[t5.distance] = kp.distance;

writetable(t5,'G:/My Drive/INVESTIGACION/POSDOC/TopoToolbox/subcuencas/rionegro_12m_knicpoint.csv','WriteMode','overwrite');
%%
% calculate theta
S = STREAMobj(FD,A>10);
mn1 = mnoptim(S,DEM,A,'mnrange',[0.1 1],'optvar','mn','crossval',false,'plot',false);
mn2 = mnoptimvar(S,DEM,A,'varfun',@robustcov);
t6 = table;
[t6.mn1] = mn1.mn;
[t6.mn2] = mn2;
writetable(t6,'G:/My Drive/INVESTIGACION/POSDOC/TopoToolbox/subcuencas/guadalupe_12m_theta.csv','WriteMode','overwrite');






