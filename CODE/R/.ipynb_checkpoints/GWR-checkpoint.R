library(spgwr)
library(sf) #simple features for spatial data
library(spdep)
library(spatialreg)
library(dplyr) #mutate


aoi = st_read("G:/My Drive/INVESTIGACION/POSDOC/Data/Vector/df_catchments_kmeans.gpkg",quiet = TRUE)
aoi$lands_dens=aoi$lands_rec/aoi$Área
aoi2 <- aoi %>% mutate_at(c('rainfallAnnual_mean','elev_mean','slope_mean','RainfallDaysmean'), ~(scale(.) %>% as.vector))
myvars <- aoi2 %>% dplyr::select(lands_dens,rainfallAnnual_mean,elev_mean,slope_mean,RainfallDaysmean)
coord=as.data.frame(st_coordinates(aoi))
coord$X
GWRbandwidth <- gwr.sel(lands_dens ~ rainfallAnnual_mean + elev_mean + slope_mean + RainfallDaysmean, 
                        data = myvars, 
                        coords=cbind(coord$X,coord$Y),
                        adapt=T)

#run the gwr model
gwr.model = gwr(lands_dens ~ rainfallAnnual_mean + elev_mean + slope_mean + RainfallDaysmean, 
                data = myvars, 
                coords=cbind(coord$X,coord$Y),
                adapt=GWRbandwidth,
                hatmatrix=TRUE,
                se.fit=TRUE)

#print the results of the model
gwr.model