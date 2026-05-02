%DEM = GRIDobj('../../Nechi/Subcuencas/ayura_12m.tif');
%FD = FLOWobj(DEM,'preprocess','carve');
A = flowacc(FD).*(FD.cellsize^2);
S = STREAMobj(FD,A>100000); 
S = klargestconncomps(S,1);
%%
imageschs(DEM,[],'colormap',[1 1 1],'colorbar',false);
hold on
plot(S)
%%
[~,mask] = maplateral(S,A,100,[]); % Map x m on each side of channel

D = DIVIDEobj(FD,S,'outlets',true,'network',false,'verbose',true);
Dc = divnet(D,FD);
D_s = sort(Dc);
D_o = divorder(D_s);
Df = removeshortdivides(D_o,FD,1000);  % Remove divides shorter than x m

% Converting DIVIDEobj into GRIDobj
[x,y]=ind2coord(Df,Df.IX);
Divides = line2GRIDobj(DEM,x,y);
Divides.Z(mask.Z==1)=0;

smooth_scale = 15; % in meters
s = smooth_scale/(sqrt(2)*pi*DEM.cellsize);

C = GRIDobj(DEM); % Prepare curvature GRIDobj
[C.Z,~,~] = conv2_rick_curv(DEM.Z,s,DEM.cellsize);
%% plot
figure
imageschs(DEM,[],'colorbar',false);
hold on
plot(Df,'color','k');
plot(S,'w');

%% plot curvature

figure
subplot(1,2,1)
imageschs(DEM,C,'caxis',[-0.1 0.1]);
 
subplot(1,2,2)
histogram(C.Z,20);
xlabel('Curvature (m^{-1})','FontSize',15);
ylabel('Counts','FontSize',15);