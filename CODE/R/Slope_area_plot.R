######SLOPE AREA PLOT

#df$drain_area <- Acumulación de flujo de la celda 
#df$Slope <- pendiente de la celda, en gradiente (si x es pendiente en grados, gradiente=((tan(x*(pi/180))))


binned <- data.frame()#matrix(ncol=3,nrow=0, dimnames=list(NULL, c('xmean', 'ymean', 'ystd'))))


#deifinición de los valores en el eje x (area)
xmin= min(df$drain_area)
xmax= max(df$drain_area)
xrange= log10(xmax)-log10(xmin) #Valores son en escala Logaritmica

#Definir el numero de intervalos en que se quiere dividir el eje x
number_of_intervals=200
interval_size=xrange/number_of_intervals

xlimits=rep(0L, number_of_intervals)
range_fun = 1:(number_of_intervals+1)

#Calculo del valor medio de cada intervalo
for (i in range_fun){
  if (i == 1){
    xlimits[1] = xmin
    } else {
      xlimits[i]= xlimits[i-1] * (10**interval_size)
    }
  }

range_fun = 1:length(xlimits)
for (j in range_fun){
  limit=xlimits[j]
  sup_limit= xlimits[j+1]
  sub <- df[df$drain_area>= limit & df$drain_area < sup_limit, ]
  xmean <- ((sup_limit+limit)/2)
  ymean= mean(sub$'Slope_Gradient')
  
  output= data.frame("xmean"= xmean, "ymean"= ymean)
  
  binned <- rbind(binned, output)}}


#Plot
ggplot(binned,aes(y=ymean,x= xmean))+ geom_point(alpha=1, size=1)+ 
labs(x = "Drainage Area (Km2)", 
   y = "Slope Gradient", 
   title = "Slope-area plots of affected basins") +
scale_y_continuous(trans='log10') +
scale_x_continuous(trans='log10')

