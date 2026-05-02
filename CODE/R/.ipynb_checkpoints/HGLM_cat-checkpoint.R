library(brms)
library(tidyverse)
library(sjstats) #for calculating intra-class correlation (ICC)
library(effects) #for plotting parameter effects
library(mosaic)  # standardizing variables
library(lme4) #The package we use, lme4. An option is to use 'nlme'.
#library(plotly)#...and this one to generate interactive plots...
library(ggplot2)#...and this one to generate interactive plots...
library(ggExtra)#to adjust plots
library(bbmle) # to get delta AICs and BICs
library(sjPlot) # more plots... especially for coeffs for multilevel models...
library(lattice) #more plots
require(MuMIn) #To get Marginal and Conditional R2
require(AICcmodavg) #To get AICc
require(pbkrtest) #... and Kenward-Roger approximations of p-values

df=read.csv("G:/My Drive/INVESTIGACION/POSDOC/GLM/df2.csv")
df$lands_dens=df$lands_rec/df$Área
df <- df[-c(4)]
df <- cbind(df, set_names(lapply(df[5:12],\(x) (x - mean(x))/sd(x)),paste0(names(df[5:12]), '_Std')))
#####################################
#functions
r2.corr.mer <- function(m) {
  lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
  summary(lmfit)$r.squared
}


######################################
fligner.test(lands_rec ~ hypso_inte, data=df) 


xyplot(slope_mean ~ hypso_inte | kmeans, data=df, type=c("p","r","smooth"))

######################################
#GLM

#Intercept model
glm_m1 <- glm(formula = lands_rec ~ 1,family = poisson,data = df)
summary(glm_m1)
coefs_glm_m1 <- data.frame(coef(summary(glm_m1)))
write.csv(coefs_glm_m1,'G:/My Drive/INVESTIGACION/POSDOC/GLM/coefs_glm_m1.csv')
AICc(glm_m1) #Corrected AIC's, more robust for small samples
logLik(glm_m1,REML=FALSE) #Higher (i.e., closer to zero) =better
deviance(glm_m1,REML=FALSE) #lower=better
r2.corr.mer(glm_m1)
1-var(residuals(glm_m1))/(var(model.response(model.frame(glm_m1))))
r.squaredGLMM(glm_m1) 
fligner.test(residuals(glm_m1) ~ df$hypso_inte) ##Fligner-Killeen's test of homogeneity of variance, one per independent variabel in each favoured model. Anything above p=0.05 is ok...

plot(fitted(glm_m1),residuals(glm_m1), family ="serif")
abline(0, 0)
lines(smooth.spline(fitted(glm_m1), residuals(glm_m1)))

#density plot of the residuals... We don't really need to check for  normality, but hey...
hist(residuals(glm_m1), prob=TRUE, ylim=c(0,0.6), family="serif")
lines(density(residuals(glm_m1)))
curve(dnorm(x, mean=mean(residuals(glm_m1)), sd=sd(residuals(glm_m1))), lty=3, col="#FF4500", add=TRUE)

#qq plot
qqnorm(residuals(glm_m1), family="serif")
qqline(residuals(glm_m1))


## Full model
glm_m2 <- glm(formula = lands_rec ~ hypso_inte_Std + X_Ksmean_Std + rainfall_cat_Std + Densidad_Std + elev_mean_Std + slope_mean_Std + rel_mean_Std,
                 family = poisson,
                 data = df)
summary(glm_m2)
coefs_glm_m2 <- data.frame(coef(summary(glm_m2)))
write.csv(coefs_glm_m2,'G:/My Drive/INVESTIGACION/POSDOC/GLM/coefs_glm_m2.csv')
AICc(glm_m2) #Corrected AIC's, more robust for small samples
logLik(glm_m2,REML=FALSE) #Higher (i.e., closer to zero) =better
deviance(glm_m2,REML=FALSE) #lower=better
r2.corr.mer(glm_m2)
1-var(residuals(glm_m2))/(var(model.response(model.frame(glm_m2))))
r.squaredGLMM(glm_m2) 


anova(glm_m1, glm_m2)

####################################
# HGLM
## Intercept model
hglm_m1 <- glmer(formula = lands_rec ~ 1 + (1|kmeans),family = poisson,data=df)
summary(hglm_m1)
coefs_hglm_m1 <- data.frame(coef(summary(hglm_m1)))
write.csv(coefs_hglm_m1,'G:/My Drive/INVESTIGACION/POSDOC/GLM/coefs_hglm_m1.csv')
AICc(hglm_m1) #Corrected AIC's, more robust for small samples
logLik(hglm_m1,REML=FALSE) #Higher (i.e., closer to zero) =better
deviance(hglm_m1,REML=FALSE) #lower=better
r2.corr.mer(hglm_m1)
1-var(residuals(hglm_m1))/(var(model.response(model.frame(hglm_m1))))
r.squaredGLMM(hglm_m1) 
performance::icc(hglm_m1)



## Full model
hglm_m2 <- glmer(formula = lands_rec ~ hypso_inte_Std + X_Ksmean_Std + rainfall_cat_Std + Densidad_Std + elev_mean_Std + slope_mean_Std + rel_mean_Std + (1|kmeans),
                    family = poisson,
                    data = df)
summary(hglm_m2)
coefs_hglm_m2 <- data.frame(coef(summary(hglm_m2)))
write.csv(coefs_hglm_m2,'G:/My Drive/INVESTIGACION/POSDOC/GLM/coefs_hglm_m2.csv')
AICc(hglm_m2) #Corrected AIC's, more robust for small samples
logLik(hglm_m2,REML=FALSE) #Higher (i.e., closer to zero) =better
deviance(hglm_m2,REML=FALSE) #lower=better
r2.corr.mer(hglm_m2)
1-var(residuals(hglm_m2))/(var(model.response(model.frame(hglm_m2))))
r.squaredGLMM(hglm_m2) 
performance::icc(hglm_m2)


## random effects
hglm_m3 <- glmer(formula = lands_rec ~ hypso_inte_Std + X_Ksmean_Std + rainfall_cat_Std + Densidad_Std + elev_mean_Std + slope_mean_Std + rel_mean_Std + (1 + hypso_inte_Std + X_Ksmean_Std + rainfall_cat_Std + Densidad_Std + elev_mean_Std + slope_mean_Std + rel_mean_Std |kmeans),
                    family = poisson,
                    data = df)
summary(hglm_m3)
coefs_hglm_m3 <- data.frame(coef(summary(hglm_m3)))
write.csv(coefs_hglm_m3,'G:/My Drive/INVESTIGACION/POSDOC/GLM/coefs_hglm_m3.csv')
AICc(hglm_m3) #Corrected AIC's, more robust for small samples
logLik(hglm_m3,REML=FALSE) #Higher (i.e., closer to zero) =better
deviance(hglm_m3,REML=FALSE) #lower=better
r2.corr.mer(hglm_m3)
1-var(residuals(hglm_m3))/(var(model.response(model.frame(hglm_m3))))
r.squaredGLMM(hglm_m3) 
performance::icc(hglm_m3)


AIC(hglm_m1,hglm_m2,hglm_m3) #lower=better. AIC good criteria for getting good prediction-models.

bbmle::AICtab(hglm_m1,hglm_m2,hglm_m3) #Call the bbmle-package first. gives delta AIC

BIC(hglm_m1,hglm_m2,hglm_m3) #lower=better. BIC usually favor simpler models. 

########################################
#plot

#To get the data right, you can use reshape2 and reshapeGUI (has to be run in R, not Rstudio).
#or that other program, begins with ex and ends with cel. I won't tell.
#e.g.
plotFULL=read.csv('G:/My Drive/INVESTIGACION/POSDOC/GLM/coefs_hglm_m3.csv') 

#Tell R the order of the variables
plotFULL$Variable <- factor(plotFULL$Variable, levels = plotFULL$Variable[order(-plotFULL$X)])

