%DEM = GRIDobj('../../Data/AOI_60m_full.tif');
%FD  = FLOWobj(DEM);
%A   = flowacc(FD);
%S   = STREAMobj(FD,A>10000);
%S   = removeedgeeffects(S,FD);
%D   = drainagebasins(FD,S);

[~,zb] = zerobaselevel(S,DEM);
%c   = zb + chitransform(S,A);
C   = mapfromnal(FD,S,c);


for r = 1:max(D)
    I = D==r;
    I = erode(I,ones(21)) & I;
    C.Z(I.Z) = nan;
end

C1 = dilate(C,ones(21));
C2 = erode(C,ones(21));

C  = C/((C1 + C2)/2);
imageschs(DEM,C,'colormap',ttscm('roma'),...
    'colorbarylabel','\chi-inferiority',...
    'ticklabels','nice');
hold on
plot(S,'k')
hold off