%% Flow direction and drainage area
DEM = GRIDobj('../../Nechi/Subcuencas/ayura_12m.tif');
DEMf = fillsinks(DEM); 
FD= FLOWobj(DEMf); 
A = flowacc(FD);  
 
%% Extract rivers
% Set the critical drainage area
res = DEM.cellsize;
A_river=1e7/(res^2); % m^2
W = A>A_river;
% Extract the wanted network
S = STREAMobj(FD,W);

%%
theta = orientation(S,'k',1000);
plotc(S,theta)
caxis([0 360])
colormap('hsv')
hx = colorbar;
hx.Label.String = 'orientation';

%1. Create a mapping structure array using STREAMobj2mapstruct
MS = STREAMobj2mapstruct(S,'attributes',{'ori' theta @mean});
%2. Export as shapefile
shapewrite(MS,'orientations.shp')

STREAMobj2kml(S,'attribute',theta,'filename','orientations.kmz','LineWidth',4);
 
%% River segments
% Extract segment and compute their orientations
mnratio=0.5;
[segment]=networksegment(S,FD,DEM,A,mnratio);
% Plot results
plotsegmentgeometry(S,segment)