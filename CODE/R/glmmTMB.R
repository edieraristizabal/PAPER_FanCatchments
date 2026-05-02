########################SPATIAL REGRESSION MODELS#####################
# datascience+ <https://datascienceplus.com/spatial-regression-in-r-part-1-spamm-vs-glmmtmb/>
#<https://backend.orbit.dtu.dk/ws/portalfiles/portal/154739064/Publishers_version.pdf>

library(tidyverse)
library(gridExtra)
library(sf) 
library(NLMR)
library(DHARMa)
library(glmmTMB)

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
col.fit1 = lm(lands_rec ~ NOMZH + elev_mean + slope_mean + RainfallDaysmean, data = aoi2)
summary(col.fit1)

# test for spatial autocorrelation
sims <- simulateResiduals(col.fit1)
testSpatialAutocorrelation(sims, aoi2$elev_mean, aoi2$lands_rec, plot = FALSE)

# fitst we need to create a numeric factor recording the coordinates of the sampled locations
aoi2$pos <- numFactor(scale(aoi2$x), scale(aoi2$y))
# then create a dummy group factor to be used as a random term
aoi2$ID <- factor(rep(1, nrow(aoi2)))

# fit the model
m_tmb <- glmmTMB(lands_rec ~ elev_mean + slope_mean + RainfallDaysmean + mat(pos + 0 | cuenca), aoi2,family = "poisson")
summary(m_tmb)

sims <- simulateResiduals(m_tmb)
plot(sims)

