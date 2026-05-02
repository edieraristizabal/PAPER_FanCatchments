################################################
##############INLA
################################################

library(sf) 
library(spdep)
library(INLA)
library(sp)
library(viridis)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(ggspatial) #scale

#Data
aoi = st_read("G:/My Drive/INVESTIGACION/POSDOC/Data/Vector/df_catchments_kmeans.gpkg",quiet = TRUE)
aoi <- aoi %>% mutate_at(c('elev_mean','slope_mean','RainfallDaysmean'), ~(scale(.) %>% as.vector))
aoi <- aoi %>% mutate(cuenca_num = case_when(cuenca == 'Atrato' ~ 1, cuenca == 'Cauca' ~ 2, cuenca == 'Magdalena' ~ 3))
aoi_sp <- as_Spatial(aoi)

ggplot() + geom_sf(data=aoi,aes(fill=lands_rec),color = "black") +
  annotation_scale(location="br",style = "ticks") +
  annotation_north_arrow(location = "tr",which_north = "true", height = unit(1, "cm"), width = unit(1, "cm"),) +
  scale_fill_gradientn(colors=c("white","orange"),name = "Landslides") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black",fill = NA,size = 1),
    axis.ticks.length=unit(-0.1, "cm"),
    axis.text.x = element_text(size = 10, margin = unit(c(t = 1, r = 0, b = 0, l = 0), "mm")),
    axis.text.y = element_text(size = 10, margin = unit(c(t = 0, r = 1, b = 0, l = 0), "mm")),
    legend.text = element_text(size=12),
    legend.title.align = 0,
    legend.position = c(0.34,0.9), 
    legend.key.size = unit(0.5, 'cm'),
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.title = element_text(size=14, vjust = .8, hjust = .5)) 

#Spatial matrix
aoi.nb <- poly2nb(aoi) #Queen matrix
# Create sparse adjacency matrix
aoi.mat <- as(nb2mat(aoi.nb, style = "B"), "Matrix")
aoi.listw = nb2listw(aoi.nb)
colnames(aoi.mat) <- rownames(aoi.mat) 
mat <- as.matrix(aoi.mat[1:dim(aoi.mat)[1], 1:dim(aoi.mat)[1]])
#Graph
nb2INLA("cl_graph",aoi.nb)
am_adj <-paste(getwd(),"G:/My Drive/INVESTIGACION/POSDOC/Figuras3/inla.graph",sep="")
H<-inla.read.graph(filename="G:/My Drive/INVESTIGACION/POSDOC/Figuras3/inla.graph")
image(inla.graph2matrix(H), xlab = "", ylab = "")
#knn5
nbs<-knearneigh(st_centroid(aoi), k = 5, longlat = T) #k=5 nearest neighbors
nbs<-knn2nb(nbs, row.names = aoi$id, sym = T) #force symmetry!!


## Standard Poisson model in glm
m1_glm = glm(lands_rec ~ RainfallDaysmean + elev_mean + slope_mean, family="poisson",offset=log(area),data = aoi)
summary(m1_glm)
aoi_sp$glm <- m1_glm$fitted.values
aoi$res=residuals(m1_glm)

ggplot() + geom_sf(data=aoi,aes(fill=res),color = "black") +
  annotation_scale(location="br",style = "ticks") +
  annotation_north_arrow(location = "tr",which_north = "true", height = unit(1, "cm"), width = unit(1, "cm"),) +
  scale_fill_gradientn(colors=c("white","orange"),name = "Residuals") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black",fill = NA,size = 1),
    axis.ticks.length=unit(-0.1, "cm"),
    axis.text.x = element_text(size = 10, margin = unit(c(t = 1, r = 0, b = 0, l = 0), "mm")),
    axis.text.y = element_text(size = 10, margin = unit(c(t = 0, r = 1, b = 0, l = 0), "mm")),
    legend.text = element_text(size=12),
    legend.title.align = 0,
    legend.position = c(0.34,0.9), 
    legend.key.size = unit(0.5, 'cm'),
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.title = element_text(size=14, vjust = .8, hjust = .5)) 

m1_over=glm(lands_rec ~ RainfallDaysmean + elev_mean + slope_mean, offset=log(area),data = aoi, family = quasipoisson)
summary(m1_over)


boxplot(res ~ cuenca, data=aoi, notch = TRUE,col=rainbow(3),xlab = "", ylab = "Residual")

## Use morans I Monte Carlo simulations
moran.mc(residuals(m1_glm),nsim = 999,listw = aoi.listw,alternative = "greater")

#Standard Poisson model in INLA
m1_inla <- inla(lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean,
           offset=log(area), data = as.data.frame(aoi),family = "poisson",
           control.predictor = list(compute = TRUE),
           control.compute = list(dic = TRUE, waic = TRUE))

summary(m1_inla)
aoi_sp$inla_m1 <- m1_inla$summary.fitted[, "mean"]

plot(m1_inla$marginals.fixed[[1]], type = "l", main = "", ylab = "", xlab = expression(beta[0]))


#basic random intercept (different intercept for each catchment)
iid_inla <- inla(lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean + 
                    f(id, model = "iid"),
                    offset=log(area), data = as.data.frame(aoi), family = "poisson", 
                    control.predictor = list(compute = TRUE),
                    control.compute = list(dic = TRUE, waic = TRUE))

summary(iid_inla)
aoi_sp$IID.inla <- iid_inla$summary.fitted[, "mean"]


#basic random intercept (cuenca)
cuenca_inla <- inla(lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean + 
           f(cuenca, model = "iid"),
           offset=log(area), data = as.data.frame(aoi), family = "poisson", 
           control.predictor = list(compute = TRUE),
           control.compute = list(dic = TRUE, waic = TRUE))

summary(cuenca_inla)
cuenca_inla$summary.random$cuenca
aoi_sp$cuenca.EFF <- cuenca_inla$summary.fitted[, "mean"]


# ICAR model
inla.icar <- inla(formula = lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean + 
                  f(cuenca, model = "iid") +
                  f(id, model = "besag", graph = aoi.mat),
                  offset=log(area), data = as.data.frame(aoi), family ="poisson",
                  control.predictor = list(compute = TRUE),
                  control.compute = list(dic = TRUE, waic = TRUE))

summary(inla.icar)
inla.icar$summary.random$cuenca
aoi_sp$ICAR <- inla.icar$summary.fitted.values[, "mean"]
 

# BYM model
inla.bym<-inla(formula = lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean +
           f(cuenca, model = "iid") +
           f(id, model = "bym",graph = aoi.mat), 
           offset=log(area), data = as.data.frame(aoi), family = "poisson",
           control.compute = list(dic = TRUE, waic=T), 
           control.predictor = list(compute = TRUE))

summary(inla.bym)
inla.bym$summary.random$cuenca
aoi$bym_id<-inla.bym$summary.random$id$mean[1:526]
aoi_sp$BYM <- inla.bym$summary.fitted.values[, "mean"]


#Leroux et al. model
ICARmatrix <- Diagonal(nrow(aoi.mat), apply(aoi.mat, 1, sum)) - aoi.mat
Cmatrix <- Diagonal(nrow(aoi), 1) -  ICARmatrix
max(eigen(Cmatrix)$values) #just to check =1

inla.ler = inla(formula = lands_rec ~ 1 + RainfallDaysmean + elev_mean + slope_mean +
                f(cuenca, elev_mean, model = "iid") +
                f(id, model = "generic1", Cmatrix = Cmatrix),
                offset=log(area), data = as.data.frame(aoi), family ="poisson",
                control.predictor = list(compute = TRUE),
                control.compute = list(dic = TRUE, waic = TRUE))

summary(inla.ler)
inla.ler$summary.random$cuenca
aoi_sp$LEROUX <- inla.ler$summary.fitted.values[, "mean"]
aoi$LEROUX <- inla.ler$summary.fitted.values[, "mean"]


###################Results###############
#spatial random effect for BYM model
ggplot() + geom_sf(data=aoi,aes(fill=bym_id),color = "black") +
  annotation_scale(location="br",style = "ticks") +
  annotation_north_arrow(location = "tr",which_north = "true", height = unit(1, "cm"), width = unit(1, "cm"),) +
  scale_fill_gradientn(colors=c("blue","red"),name = "Spatial random effect") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black",fill = NA,size = 1),
    axis.ticks.length=unit(-0.1, "cm"),
    axis.text.x = element_text(size = 10, margin = unit(c(t = 1, r = 0, b = 0, l = 0), "mm")),
    axis.text.y = element_text(size = 10, margin = unit(c(t = 0, r = 1, b = 0, l = 0), "mm")),
    legend.text = element_text(size=12),
    legend.title.align = 0,
    legend.position = c(0.34,0.9), 
    legend.key.size = unit(0.5, 'cm'),
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.title = element_text(size=14, vjust = .8, hjust = .5)) 

#
ggplot(data = aoi)+geom_histogram(aes(x =lands_rec,y=..density..))


spplot(aoi_sp, c("lands_rec", "ICAR", "BYM", "LEROUX"),col.regions = rev(magma(16)))
spplot(aoi_sp, c("IID.inla", "ICAR", "BYM", "LEROUX"),col.regions = rev(magma(16)))
spplot(aoi_sp, c("lands_rec", "glm","inla_m1","IID.inla"),col.regions = rev(magma(16)))


