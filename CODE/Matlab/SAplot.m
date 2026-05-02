DEM = GRIDobj('G:/My Drive/INVESTIGACION/POSDOC/Data/Raster/subcuencas/dem_V5_CUENCA292_292.tif');
FD = FLOWobj(DEM,'preprocess','carve');
A = flowacc(FD);
S = STREAMobj(FD,A>100);
DEM = imposemin(S,DEM);
%S = removeshortstreams(S,100);
%S   = trunk(S);

a = getnal(S,A)*(DEM.cellsize^2);
g = gradient(S,DEM);

label = labelreach(S,'seglength',1000);

gg = accumarray(label,g,[],@mean);
ggs = accumarray(label,g,[],@(x) std(x)/sqrt(numel(x)));

ag = accumarray(label,a,[],@mean);
ags = accumarray(label,a,[],@(x) std(x)/sqrt(numel(x)));

% plot
plot([ag ag]',[gg+ggs max(gg-ggs,0.0001)]','color',[.7 .7 .7]);
hold on
plot([ag-ags ag+ags]',[gg gg]','color',[.7 .7 .7]);
plot(ag,gg,'o','MarkerFaceColor',[.7 .7 .7])
hold off
set(gca,'Xscale','log','Yscale','log');
xlabel('Area [m^2]')
ylabel('Gradient [m m^{-1}]')
ylim([1e-3 0.6])