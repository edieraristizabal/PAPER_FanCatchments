########################SPATIAL REGRESSION MODELS#####################
# datascience+ <https://datascienceplus.com/spatial-regression-in-r-part-1-spamm-vs-glmmtmb/>

library(tidyverse)
library(gridExtra)
library(sf) 
library(NLMR)
library(DHARMa)
library(spaMM)

aoi = st_read("G:/My Drive/INVESTIGACION/POSDOC/Data/Vector/df_catchments_kmeans.gpkg",quiet = TRUE)
aoi2 <- aoi %>% mutate_at(c('elev_mean','slope_mean','RainfallDaysmean'), ~(scale(.) %>% as.vector))
aoi2$lands_dens=aoi$lands_rec/aoi$Área

## Store the geometry of the polygons
aoi.geom = st_geometry(aoi)

## Store the centroids of the polygons
aoi.coords = st_centroid(aoi.geom)

#to get coordinates
aoi_xy=st_coordinates(aoi.coords)
aoi2$x=aoi_xy[1:533,1]
aoi2$y=aoi_xy[1:533,2]

## Basic linear model
col.fit1 = lm(lands_rec ~ elev_mean + slope_mean + RainfallDaysmean, data = aoi2)
summary(col.fit1)

# test for spatial autocorrelation
sims <- simulateResiduals(col.fit1)
testSpatialAutocorrelation(sims, aoi2$elev_mean, aoi2$lands_dens, plot = FALSE)

# fit the spatial model
m_spamm <- fitme(lands_rec ~ NOMZH + elev_mean + slope_mean + RainfallDaysmean + Matern(1 | x + y), data = aoi2, family = "poisson")
summary(m_spamm)

dd <- dist(cbind(as.vector(aoi2$x),as.vector(aoi2$y)))
mm <- MaternCorr(dd, nu = 0.270, rho = 0.00002)
plot(as.numeric(dd), as.numeric(mm), xlab = "Distance between pairs of location [in m]", ylab = "Estimated correlation")


sims <- simulateResiduals(m_spamm)
plot(sims)  
