
################ In this R file, we demonstrate the use of the developed methodology and reproduce
################ the estimation results of the empirical land price study.
################ Note that because the model is implemented using the Bayesian MCMC approach, the estimation
################ results of the land price study might be slightly different from Table 2 in the manuscript

rm(list=ls())

### set working directory for R
setwd("M:/Practical 2/Practical 2/")
### load necessary R packages to use the function that implement the developed model
library(MCMCpack) 
library(mvtnorm) 
library(spdep)   
library(spatialprobit) 
library(Matrix)
library(foreign)
library(rgdal)
library(RColorBrewer)

### Import the land price data set used in the paper
data <- read.csv("landpricedata.csv",header=TRUE)

### Order the data according to the variable, streetid, that groups land parcels into districts

data <- data[order(data$streetid),] 
MM <- as.data.frame(table(data$streetid))
Utotal <- dim(MM)[1]
Unum <- MM[,2]
Uid <- rep(c(1:Utotal),Unum)

### Extract the dependent variable y and independent variables X
# dependent variable y
y <- data$y
# independent variables
varnames <- c("Intercept","Logarea","LogDcbd","LogDsubway","LogDele","LogDpark","LogDriver",
              "Popden","Buildings1949","Crimerate","Year04","Year05","Year06","Year07",
              "Year08","Year09")
X <- data[varnames]
X <- as.matrix(X)
n <- dim(X)[1]
p <- dim(X)[2]

### Import other arguments in the function that implement the developed model
# The land parcel level (lower level) spatial weights matrix W
load("W.RData")
# The district level (higher level) spatial weights matrix M
load("M.RData")
# The random effect design matrix Z
load("Z.RData")

# the Log-determinants of two Jacobian terms 
tmp <- sar_lndet(ldetflag=2,W,rmin,rmax)
detval <- tmp$detval
tmp <- sar_lndet(ldetflag=2,M,rmin,rmax)
detvalM <- tmp$detval

####### Running the function names HSAR.R
library(distr) 
source("HSAR.R")
res1 <- HSAR(y=y,X=X,W=W,M=M,Z=Z)

####### Now the estimation results, see Table 2 in the manuscript
res1
# The regression coefficients
dd <- data.frame(res1$Mbetas,res1$SDbetas,t=res1$Mbetas/res1$SDbetas)
row.names(dd) <- varnames
dd

# rho, land parcel level spatial autoregressive parameter
data.frame(res1$Mrho,res1$SDrho,t=res1$Mrho/res1$SDrho)

# lambda, the district level spatial autoregressive parameter
data.frame(res1$Mlambda,res1$SDlambda,t=res1$Mlambda/res1$SDlambda)

# land parcel level variance sigma2e
data.frame(res1$Msigma2e,res1$SDsigma2e,t=res1$Msigma2e/res1$SDsigma2e)

# district parcel level variance sigma2u
data.frame(res1$Msigma2u,res1$SDsigma2u,t=res1$Msigma2u/res1$SDsigma2u)
