DEM = GRIDobj('../../Data/AOI_60m_full.tif');
[x,y,p] = prominence(DEM,500);
imageschs(DEM,[],'colormap',[1 1 1],'colorbar',false);
hold on
h = bubblechart(x,y,p);
h.MarkerEdgeColor = 'k';
bubblelim([200 2000])
bubblesize([2 20])
bubblelegend('location','southeast')