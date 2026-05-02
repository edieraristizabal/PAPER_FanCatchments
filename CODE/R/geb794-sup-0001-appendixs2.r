# R code supplement to:
# Swanson, A, Dobrowski, S, Finley, A., Thorne, J. & Schwartz, M. (2012) Spatially explicit 
# methods capture prediction uncertainty in species distribution model forecasts through time.
# Global Ecology and Biogeography, XX, XXX-XXX.
#<https://onlinelibrary.wiley.com/doi/10.1111/j.1466-8238.2012.00794.x>
# This script shows how to generate artificial data, fit a spatial process GLMM, run diagnostics, 
# assess prediction accuracy and generate gridded predictions.
#
# Required libraries (available through CRAN) include spBayes, geoR, PresenceAbsence, pgirmess, 
#  and fields.
# 
# Author: Alan Swanson (alan.swanson@umontana.edu)
###############################################################################


#############################################################################################
# generate synthetic data.
#############################################################################################
library(spBayes)
library(geoR)
logit <- function(x) 1/(1+exp(-x))
N <- 50 # specify a 50 by 50 grid of data
n <- 250 # randomly sample 250 points from this grid.  
true_coef <- c(-1,3) # true coefficient values
nsamp <- rep(10,n) # number of samples per observation

# generate a random intercept.
set.seed(123456)
phi <- 3/25 # for exponential covariance, this is equivelent to an effective range of ~25 units.
sigma.sq <- 1 # sd=1
w <- grf(N^2,grid="reg",xlim=c(1,N),ylim=c(1,N),cov.model="exponential",cov.pars=c(sigma.sq,1/phi),messages=F)

# generate a random covariate.
phi_x <- 3/25
sigma_x <- 1
x1 <- grf(N^2,grid="reg",xlim=c(1,N),ylim=c(1,N),cov.model="exponential",cov.pars=c(sigma_x,1/phi_x),messages=F)

# combine to get true probability of occurrence.
p <- logit(cbind(rep(1,N^2),x1$data) %*% matrix(true_coef,ncol=1)+w$data)
			
# sample a random subset.
ss <- sample(N^2,n) # random subset
ypres <- rbinom(n,nsamp,p[ss]) # generate response
yabs <- nsamp-ypres  # absences
x1samp <- x1$data[ss] # covariate values at sampled locations

# fit a GLM to get starting values
m0 <- glm(cbind(ypres,yabs)~x1samp,family="binomial") 
	
# spatial process GLMM starting/tuning parameters #
beta.starting <- coefficients(m0)
beta.tuning <- t(chol(vcov(m0)))
n.batch <- 500
batch.length <- 50
n.samples <- n.batch*batch.length

# knot coords
k.coords <- kmeans(w$coords[ss,],50)$centers
k.dist <- iDist(k.coords)
min.dist <- apply(k.dist,1,function(x) min(x[x>0]))
summary(min.dist)  # mean nearest-neighbor dist between knots is approx. 5

# spatial process GLMM starting/tuning parameters #
beta.starting <- coefficients(m0)
beta.tuning <- t(chol(vcov(m0)))
n.batch <- 500
batch.length <- 50
n.samples <- n.batch*batch.length
phi.lims <- 3/c(50,5) # limit domain of phi based on inter-knot distances.
phi.strt <- 3/20
sigma.strt <- 0.5

#############################################################################################
# fit three independent chains of spatial GLMM model.
# this stage will take more than an hour on most computers.
#############################################################################################
m1a <- spGLM(ypres~x1samp, family="binomial", coords=w$coords[ss,],weights=nsamp, 
	knots=k.coords,  # by giving knot coordinates we specify a predictive process model.
	starting=list("beta"=beta.starting, "phi"=phi.strt,"sigma.sq"=sigma.strt, "w"=0,"e"=0),
	tuning=list("beta"=beta.tuning, "phi"=0.5, "sigma.sq"=0.1, "w"=0.01,"e"=0),
	priors=list("beta.Normal"=list(c(0,0),c(10,10)),"phi.Unif"=phi.lims,"sigma.sq.IG"=c(2, 1)),
	amcmc=list("n.batch"=n.batch,"batch.length"=batch.length,"accept.rate"=0.43),
	cov.model="exponential",
	modified.pp=T,  # specifying a modified predictive process 
	n.samples=n.samples, sub.samples=c(1,n.samples,10),
	verbose=TRUE, n.report=15)

m1b <- spGLM(ypres~x1samp, family="binomial", coords=w$coords[ss,], 
	weights=nsamp,knots=k.coords,
	starting=list("beta"=beta.starting, "phi"=phi.strt,"sigma.sq"=sigma.strt, "w"=0,"e"=0),
	tuning=list("beta"=beta.tuning, "phi"=0.5, "sigma.sq"=0.1, "w"=0.01,"e"=0),
	priors=list("beta.Normal"=list(c(0,0),c(10,10)),"phi.Unif"=phi.lims,"sigma.sq.IG"=c(2, 1)),
	amcmc=list("n.batch"=n.batch,"batch.length"=batch.length,"accept.rate"=0.43),
	cov.model="exponential",
	modified.pp=T,  
	n.samples=n.samples, sub.samples=c(1,n.samples,10),
	verbose=TRUE, n.report=15)

m1c <- spGLM(ypres~x1samp, family="binomial", coords=w$coords[ss,], 
	weights=nsamp,knots=k.coords,
	starting=list("beta"=beta.starting, "phi"=phi.strt,"sigma.sq"=sigma.strt, "w"=0,"e"=0),
	tuning=list("beta"=beta.tuning, "phi"=0.5, "sigma.sq"=0.1, "w"=0.01,"e"=0),
	priors=list("beta.Normal"=list(c(0,0),c(10,10)),"phi.Unif"=phi.lims,"sigma.sq.IG"=c(2, 1)),
	amcmc=list("n.batch"=n.batch,"batch.length"=batch.length,"accept.rate"=0.43),
	cov.model="exponential",
	modified.pp=T,  
	n.samples=n.samples, sub.samples=c(1,n.samples,10),
	verbose=TRUE, n.report=15)


#############################################################################################
# examine output from spatial GLMM models #
#############################################################################################

# Consolodate posteriors to perform diagnostics of spGLM fit.
n.samp <- nrow(m1a$p.samples)
strt <- ceiling(n.samp/2)
posteriors <- as.mcmc.list(list(m1a$p.samples,m1b$p.samples,m1c$p.samples))

# The following step converts the spatial decay parameter, phi, to effective range.
for(i in 1:3) posteriors[[i]][,"phi"]<-3/posteriors[[i]][,"phi"]  
plot(window(posteriors,start=strt))

# Compute Gelman diagnostics to assess convergence.  These compare within-chain to 
# between-chain variability. Values near 1 suggest full convergence.
gelman.diag(window(posteriors,start=strt,end=n.samp,thin=1))

# Extract median predictions at sample points.
np<-m1a$p
pooled.samples <- rbind(m1a$p.samples[strt:n.samp,1:np],m1a$p.samples[strt:n.samp,1:np],m1a$p.samples[strt:n.samp,1:np])
pooled.ranef <- cbind(m1a$sp.effects[,strt:n.samp],m1b$sp.effects[,strt:n.samp],m1c$sp.effects[,strt:n.samp])
Xb <- cbind(1,x1samp) %*% t(pooled.samples)
yhat <- logit(Xb + pooled.ranef)
yhat <- apply(yhat,1,median)

# Calculate AUC #
library(PresenceAbsence)
y.pa <- as.numeric(ypres>0) # reduce response to presence-absence.
auc(cbind(ss,y.pa,yhat))  # for GLMM
auc(cbind(ss,y.pa,fitted(m0,type="response"))) # for GLM

# Compare parameter estimates and intervals
alpha <- 0.05
coefs <- as.matrix(summary(m0)$coefficients[,1:2])
ints <- coefs%*%rbind(c(1,1),qnorm(c(alpha/2,1-alpha/2)))
t(cbind(lwr=ints[,1],est=coefs[,1],upr=ints[,2])) 
apply(pooled.samples,2,function(x) quantile(x,c(alpha/2,0.5,1-alpha/2)))
# Note that GLM intervals do not include the true values of the slope and intercept (-1 and 3 respectively)

# Calculate GLMM deviance residuals.
sgn <- ifelse(ypres-yhat*10>=0,1,-1)
z <- sgn*sqrt(2*(ypres*log(ypres/(nsamp*yhat))+(nsamp-ypres)*log((nsamp-ypres)/(nsamp*(1-yhat)))))
z[ypres==0]<- -1*sqrt(2*nsamp*abs(log(1-yhat)))[ypres==0]
z[ypres==nsamp]<- sqrt(2*nsamp*abs(log(yhat)))[ypres==nsamp]

# Calculate Moran's I statistic on deviance residuals by distance class.
library(pgirmess)
par(mfrow=c(1,2))
morans.glmm <- correlog(w$coords[ss,],z,nbclass=10)
plot(morans.glmm,main="GLMM")
abline(h=0,col=2,lty=2)
print(morans.glmm)

morans.glm <- correlog(w$coords[ss,],residuals(m0,type="deviance"),nbclass=10)
plot(morans.glm,main="GLM")
abline(h=0,col=2,lty=2)
print(morans.glm)
# Note that GLM residuals show significant autocorrelation in the first two distance classes,
# but GLMM residuals show no significant autocorrelation.

# Make predictions across the entire grid.  Note that this step requires ~2Gb of RAM to be 
# available in R.  The function memory.limit() may be needed to make more RAM available. 
pred.coords <- expand.grid(1:N,1:N)
memory.limit(6000)
preds1 <- spPredict(m1a,pred.coords=pred.coords,pred.covars=cbind(1,x1$data),start=strt,thin=10)
preds2 <- spPredict(m1b,pred.coords=pred.coords,pred.covars=cbind(1,x1$data),start=strt,thin=10)
preds3 <- spPredict(m1c,pred.coords=pred.coords,pred.covars=cbind(1,x1$data),start=strt,thin=10)

# Pool predictions across chains, find median values, and put in matrix form.
ypred <- t(matrix(apply(cbind(preds1$y.pred,preds2$y.pred,preds3$y.pred),1,median,na.rm=T),ncol=N,byrow=T))
wpred <- t(matrix(apply(cbind(preds1$w.pred,preds2$w.pred,preds3$w.pred),1,median,na.rm=T),ncol=N,byrow=T))

# Make image of fits and true values.
library(fields)  # the image() function in the base graphics package will work if fields is unavailable.
par(mfrow=c(2,2))
image.plot(x=1:N,y=1:N,z=t(matrix(p,nrow=N,byrow=T)),xlab="x",ylab="y",main="true response",zlim=c(0,1))
image.plot(x=1:N,y=1:N,z=ypred,xlab="x",ylab="y",main=expression(paste("GLMM ",hat(y))),zlim=c(0,1))
image.plot(x=1:N,y=1:N,z=t(matrix(w$data,nrow=N,byrow=T)),main="true w",xlab="x",ylab="y",zlim=c(-3.5,3.5))
image.plot(x=1:N,y=1:N,z=wpred,main=expression(paste("GLMM ",hat(w))),xlab="x",ylab="y",zlim=c(-3.5,3.5))


