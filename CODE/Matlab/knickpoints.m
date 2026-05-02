%% prepare DEM
DEM = GRIDobj('G:/My Drive/INVESTIGACION/POSDOC/Data/Raster/subcuencas/dem_V5_CUENCA219_219.tif');
FD = FLOWobj(DEM,'preprocess','carve');
A   = flowacc(FD);
S   = STREAMobj(FD,A>100);
Sp = trunk(S); %para obtener solo el drenaje principal

%%
imageschs(DEM,[],'colormap',[1 1 1],'colorbar',false);
hold on
plot(S)

%% Estimate tolerance value
zmax = cummaxupstream(S,DEM);
zmin = imposemin(S,getnal(S,DEM));
tol=max([zmax-zmin])
plotdzshaded(S,[zmin zmax])
hold on
plotdz(S,DEM,'color','k')

%% Calculate knickpoints
[~,kp] = knickpointfinder(S,DEM,'tol',20,'split',false,'plot',true);

hold on
scatter(kp.distance,kp.z,kp.dz,'sk','MarkerFaceColor','r')
hold off

%% plot map with bubbles
P = PPS(S,'PP',kp.IXgrid,'z',DEM);
imageschs(DEM,[],'colormap',[1 1 1])
hold on
plot(S)
hold on
plotpoints(P,'sizedata',kp.dz,'colordata',DEM)
axis image
hold off
xlabel('Easting [m]');
ylabel('Northing [m]');
bubblelegend('Location','Northwest')

%% plot with a DEM
imageschs(DEM)
hold on
plot(S,'k');
plot(kp.x,kp.y,'ko','MarkerFaceColor','w')
hold off

%% Plot a longitudinal river profile
plotdz(S,DEM);
hold on
plot(kp.distance,kp.z,'ko','MarkerFaceColor','w')
xlabel('\chi') 
hold off
%% plot a longitudinal river profile in chi space
A = flowacc(FD);
c = chitransform(S,A,'mn',0.45);
plotdz(S,DEM,'distance',c);
hold on
[~,locb] = ismember(kp.IXgrid,S.IXgrid);
ckp = c(locb);
plot(ckp,kp.z,'ko','MarkerFaceColor','w')
hold off
xlabel('\chi [m]')

%% Generating SHP Kinckpoints
nr_Kp = numel(kp.n); kinck = [];
for i=1:kp.n
kinck(i).Geometry = deal('Point');
kinck(i).id = i;
kinck(i).X = kp.x(i);
kinck(i).Y = kp.y(i);
kinck(i).distan = kp.distance(i);
kinck(i).z = kp.z(i);
kinck(i).IXgrid = kp.IXgrid(i);
kinck(i).order = kp.order(i);
kinck(i).dz = kp.dz(i);
end
fprintf('ready\n')

% Export knickpoints to shapefiles
%shapewrite(kinck,'Knicpoints2_100m'); 
