######### A formal/simple function for a HSAR model with higher level dependence 


### Key model inputs/parameters
#@param y, the dependent variable
#@param X, the independent variables
#@param W, the lower level spatial weight matrix, need to be row-normalized
#@param M, the higher level spatial weight matrix, need to be row-normalized
#@param detvalW, the pre-calculated LOG-determinant for W, i.e., LOG|A| in the HSAR paper,see Barry and Pace (1999)
#@param detvalM, the pre-calculated LOG-determinant for M, i.e., LOG|B| in the HSAR paper,see Barry and Pace (1999)
#@param Z, the design matrix for higher level random effects



#########################################################
#### function to implement a HSAR model #################
#########################################################

HSAR <- function(y,X,W,M,Z) {

### check for if spatial weight matrices W and M are valid

if (!inherits(W,"sparseMatrix") || any(diag(W) != 0)) {
	stop("W should be of sparse matrix form and diagonal elements should be 0s")	
}

if (!inherits(M,"sparseMatrix") || any(diag(M) != 0)) {
	stop("M should be of sparse matrix form and diagonal elements should be 0s")	
}

### change the dense random effect design matrix to a sparse one (optional)
Z <- as(Z,"dgCMatrix") 

start.time <- Sys.time()

### Some key functions in the estimation

### faster update of matrix A = I - rho*W for new values of rho
# @param A template matrix of (I - rho*W)
# @param ind indices to be replaced
# @param W spatial weight matrix with sparse matrix form
# @return (I - rho*W)

update_A <- function(A,ind,rho,W) {
	A@x[ind] <- (-rho*W)@x
	return(A)
}
###

### faster update of matrix B = I - lambda*M for new values of rho
# @param B template matrix of (I - lambda*M)
# @param ind indices to be replaced
# @param M spatial weight matrix with sparse matrix form
# @return (I - lambda*M)

update_B <- function(B,ind,lambda,M) {
	B@x[ind] <- (-lambda*M)@x
	return(B)
}
###


#########################################################
####Supporting functions for drawing rho    #############
#########################################################

### draw rho-values using griddy Gibbs and inversion, see Smith and LeSage (2004)###

draw_rho <- function(detval,e0e0,eded,eueu,e0ed,e0eu,edeu,sig) {
	rho_exist <- detval[,1]
	nrho <- length(rho_exist)
	log_detrho <- detval[,2]
	iota <- rep(1,times=nrho)
	#####Calculate Log-S(e(rho*)) given both rho* and u*
	S_rho <- e0e0*iota + rho_exist^2*eded + eueu - 2*rho_exist*e0ed - 2*e0eu + 2*rho_exist*edeu 

	##### Calculate the Log-density
	log_den <- log_detrho - S_rho/(2*sig) 
	adj <- max(log_den)
	log_den <- log_den-adj

	##### the density
	den <- exp(log_den)
	##### the interval separating rho is h=0.001
	h <- 0.001

	##### Integration to calculate the normalized constant
	##### using the  trapezoid rule

	ISUM <- h*(den[1]/2 + sum(den[2:1890]) + den[1891]/2)
	norm_den <- den/ISUM
	##### cumulative density
	cumu_den <- cumsum(norm_den)

	##### Inverse sampling
	rnd <- runif(1,0,1)*sum(norm_den)
	ind <- which(cumu_den <= rnd)
	idraw <- max(ind)
	if(idraw > 0 & idraw < nrho) rho_draw <- rho_exist[idraw]
	return(rho_draw)
}

#########################################################
####Supporting functions for drawing lambda #############
#########################################################

### draw lambda-values using griddy Gibbs and inversion, see Smith and LeSage (2004)###

draw_lambda <- function(detval,uu,uMu,uMMu,sig) {
	lambda_exist <- detval[,1]
	nlambda <- length(lambda_exist)
	log_detlambda <- detval[,2]
	iota <- rep(1,times=nlambda)

	
	S_lambda <- uu*iota - 2*lambda_exist*uMu + lambda_exist^2*uMMu

	##### Calculate the Log-density
	log_den <- log_detlambda - S_lambda/(2*sig) 
	adj <- max(log_den)
	log_den <- log_den-adj

	##### the density
	den <- exp(log_den)
	##### the interval separating lambda is h=0.001
	h <- 0.001

	##### Integration to calculate the normalized constant
	##### using the  trapezoid rule

	ISUM <- h*(den[1]/2 + sum(den[2:1890]) + den[1891]/2)
	norm_den <- den/ISUM
	##### cumulative density
	cumu_den <- cumsum(norm_den)

	##### Inverse sampling
	rnd <- runif(1,0,1)*sum(norm_den)
	ind <- which(cumu_den <= rnd)
	idraw <- max(ind)
	if(idraw > 0 & idraw < nlambda) lambda_draw <- lambda_exist[idraw]

	return(lambda_draw)
}


#######################################################
####Starting the MCMC SET UP              #############
#######################################################

#### Uniform distribution for the spatial auto-regressive parameters rho, and lambda
 

n <- dim(X)[1]
p <- dim(X)[2]


################ Prior distribution specifications

### For Betas

M0 <- rep(0,times=p)
T0 <- diag(100,p)

### For sigma2e and sigma2u 
### completely non-informative priors

c0=d0=a0=b0=0.01

################ Store MCMC results
Nsim <- 10000
burnin <- 5000
Betas <- matrix(0,nrow=Nsim,ncol=p)
Us <- matrix(0,nrow=Nsim,ncol=Utotal)
sigma2e <- rep(0,times=Nsim)
sigma2u <- rep(0,times=Nsim)
rho <- rep(0,times=Nsim)
lambda <- rep(0,times=Nsim)

#### initial values for model parameters (better initial values will be that from a classic multilevel model)
sigma2e[1] <- 2
sigma2u[1] <- 2
rho[1] <- 0.5
lambda[1] <- 0.5

################ Fixed posterior hyper-parameters, 
################ the shape parameter in the Inverse Gamma distribution

ce <- n/2 + c0
au <- Utotal/2 + a0


################ Fixed matrix manipulations during the MCMC loops

XTX <- crossprod(X)
invT0 <- solve(T0)
T0M0 <- invT0%*%M0
tX <- t(X)
tZ <- t(Z)
################# some fixed values when updating rho

beta0 <- solve(X,y)
e0 <- y-X%*%beta0
e0e0 <- crossprod(e0)

Wy <- as.numeric(W%*%y)
betad <- solve(X,Wy)
ed <- Wy-X%*%betad
eded <- crossprod(ed)

e0ed <- crossprod(e0,ed)


################# initialise A

if (class(W) == "dgCMatrix") {
	I <- sparseMatrix(i=1:n,j=1:n,x=Inf)
	A <- I - rho[1]*W
	ind <- which(is.infinite(A@x))
	ind2 <- which(!is.infinite(A@x))
	A@x[ind] <- 1
}

################# initialise B

if (class(M) == "dgCMatrix") {
	I <- sparseMatrix(i=1:Utotal,j=1:Utotal,x=Inf)
	B <- I - lambda[1]*M
	indB <- which(is.infinite(B@x))
	ind2B <- which(!is.infinite(B@x))
	B@x[indB] <- 1
}

##########################################
####MCMC updating               ##########
##########################################

for(i in 2:Nsim) {

################################# Gibbs sampler for updating Betas

VV <- XTX/sigma2e[i-1] + invT0
#### We use Cholesky decomption to inverse the covariance matrix
vBetas <- chol2inv(chol(VV))
if(inherits(vBetas,"try-error")) vBetas <- solve(VV)

### Define A=I-rho*W
A <- update_A(A,ind=ind2,rho=rho[i-1],W=W)
Ay <- as.numeric(A%*%y)
ZUs <- as.numeric(Z%*%Us[i-1,])
mBetas <- vBetas%*%(tX%*%(Ay-ZUs)/sigma2e[i-1]+T0M0)

#### When the number of independent variables is large, drawing from multivariate 
#### distribution could be time-comsuming
cholV <- t(chol(vBetas))
#### draw form p-dimensitional independent norm distribution
betas <- rnorm(p)
Betas[i,] <- betas <- mBetas + cholV%*%betas


################################# Gibbs sampler for updating U. Now, us are spatially dependent.

### Define and update B=I-lambda*M

B <- update_B(B,ind=ind2B,lambda=lambda[i-1],M=M)

vU <- as.matrix(tZ%*%Z/sigma2e[i-1] + t(B)%*%B/sigma2u[i-1])
vU <- chol2inv(chol(vU))
if(inherits(vU,"try-error")) vU <- solve(vU)

Xb <- X%*%betas
mU <- vU%*%(as.numeric(tZ%*%(Ay - Xb))/sigma2e[i-1])

#### When the number of higher level units is large, drawing from multivariate 
#### distribution could be time-comsuming
cholV <- t(chol(vU))

#### draw form J-dimensitional independent norm distribution
us <- rnorm(Utotal)
Us[i,] <- us <- mU + cholV%*%us


################################# Gibbs sampler for updating sigma2e
Zu <- rep(us,Unum) #### Z%*%us
e <- Ay-Zu-Xb
de <- 0.5*crossprod(e)+d0
sigma2e[i] <- rinvgamma(1,shape=ce,scale=de)


################################# Gibbs sampler for updating sigma2u
Bus <- as.numeric(B%*%us)
bu <- crossprod(Bus)/2 + b0
sigma2u[i] <- rinvgamma(1,shape=au,scale=bu)


################################# Giddy Gibbs integration and inverse sampling for rho

betau <- solve(X,Zu)
eu <- Zu-X%*%betau
eueu <- crossprod(eu)
e0eu <- crossprod(e0,eu)
edeu <- crossprod(ed,eu)

#### Using the draw_rho function to update rho
rho[i] <- draw_rho(detval=detval,e0e0=e0e0,eded=eded,eueu=eueu,e0ed=e0ed,e0eu=e0eu,
		edeu=edeu,sig=sigma2e[i])


################################# Giddy Gibbs integration and inverse sampling for lambda

uu <- crossprod(us)
uMu <- as.numeric(t(us) %*% M %*% us)
Mu <- as.numeric(M %*% us)
uMMu <- crossprod(Mu)

lambda[i] <- draw_lambda(detval=detvalM,uu=uu,uMu=uMu,uMMu=uMMu,sig=sigma2u[i])


if(i%%50==0) cat(".")

}

end.time <- Sys.time()
######### finish the MCMC loops

######### store the results

results <- NULL
results$time <- end.time - start.time
results$Mbetas <- apply(Betas[burnin:Nsim,],2,mean)
results$SDbetas <- apply(Betas[burnin:Nsim,],2,sd)
results$Mrho <- mean(rho[burnin:Nsim])
results$SDrho <- sd(rho[burnin:Nsim])
results$Mlambda <- mean(lambda[burnin:Nsim])
results$SDlambda <- sd(lambda[burnin:Nsim])
results$Msigma2e <- mean(sigma2e[burnin:Nsim])
results$SDsigma2e <- sd(sigma2e[burnin:Nsim])
results$Msigma2u <- mean(sigma2u[burnin:Nsim])
results$SDsigma2u <- sd(sigma2u[burnin:Nsim])
results$Mus <- apply(Us[burnin:Nsim,],2,mean)
results$SDus <- apply(Us[burnin:Nsim,],2,sd)

return(results)

}



