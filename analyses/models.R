#Graph Learning models
#Charley Wu

packages <- c("MASS", 'matrixcalc', 'Matrix', 'igraph', 'expm')
lapply(packages, require, character.only = TRUE)
source('utilities.R')

##############################################################################################################
#Diffusion Kernel
##############################################################################################################

diffusionKernel<- function(network, alpha = 1, precision = 10, normalizedLaplacian = F){
  #1. Construct Graph Primatives
  D <- diag(strength(network)) #degree matrix
  #(weighted) adjacency matrix
  if (is_weighted(network)){# Weighted graph
    W <- as_adjacency_matrix(network, attr="weight", sparse=F)
  }else{
    W <- as_adjacency_matrix(network, sparse=F)
  }
  #2. construct laplacian
  L <- D - W
  if (normalizedLaplacian==T){
    L <- sqrt(solve(D))%*%L%*%sqrt(solve(D)) #if degree normalized Laplacian
  }
  #3. Eigenvalue decomposition
  eig <- eigen(-L)
  #4. Matrix Exponential of Laplacian matrix
  K <- eig$vectors %*% diag(exp(alpha * eig$values)) %*% t(eig$vectors) #exp(alpha * L)
  K <- round(K, precision) #necessary since loss of precision can make other functions not recognize K as positive semi-definite 
  return(K)
}

class(diffusionKernel)<- c(class(diffusionKernel), 'DF')

##############################################################################################################
#Graph regularization operators: Smola & Kondor, (2003)
##############################################################################################################

#A variant of this is used in Tenenbaum & Griffiths (2009)
regularizedLaplacian <- function(network, sigma, precision = 10, normalizedLaplacian=F){
  #1. Construct Graph Primatives
  D <- diag(strength(network)) #degree matrix
  #(weighted) adjacency matrix
  if (is_weighted(network)){# Weighted graph
    W <- as_adjacency_matrix(network, attr="weight", sparse=F)
  }else{
    W <- as_adjacency_matrix(network, sparse=F)
  }
  #2. construct laplacian
  L <- D - W
  if (normalizedLaplacian==T){
    L <- sqrt(solve(D))%*%L%*%sqrt(solve(D)) #if degree normalized Laplacian
  }
  #Compute kernel
  K <- solve(diag(length(V(network))) + sigma^2 * L)
  K <- round(K, precision) #necessary since loss of precision can make other functions not recognize K as positive semi-definite 
  return(K)
}

#The diffusion kernel is the continuous time limit of this kernel, which uses random walks of fixed step lengths
randomWalkKernel <- function(network, p, A = 2,  precision = 10, normalizedLaplacian=T){
  #1. Construct Graph Primatives
  D <- diag(strength(network)) #degree matrix
  #(weighted) adjacency matrix
  if (is_weighted(network)){# Weighted graph
    W <- as_adjacency_matrix(network, attr="weight", sparse=F)
  }else{
    W <- as_adjacency_matrix(network, sparse=F)
  }
  #2. construct laplacian
  L <- D - W
  if (normalizedLaplacian==T){
    L <- sqrt(solve(D))%*%L%*%sqrt(solve(D)) #if degree normalized Laplacian
  }
  #Compute kernel
  K <-  matrix.power(A*diag(length(V(network))) - L, p)
  K <- round(K, precision) #necessary since loss of precision can make other functions not recognize K as positive semi-definite 
  return(K)
}


##############################################################################################################
#Graph predictions: Experiment 1
##############################################################################################################

#GP prediction
gpPrediction <- function(pars, kernel, dataset, normalizedLaplacian = F){
  #extract parameters
  if (kernel == 'DF'){
    alpha = pars[1]
  }else if (kernel == 'RL'){
    sigma = pars[1]
  }else if (kernel == 'RW'){
    p = pars[1]
  }
  prediction <-rep(NA, nrow(dataset)) 
  uncertainty <-rep(NA, nrow(dataset)) 
  for (i in 1:nrow(dataset)){
    row <- dataset[i,] #subset data
    g <- graphList[[row$graphOrder+1]] #acquire the graph for this trial; +1 converts from base0 to base1
    obs <- unlist(row[,c('observation1', 'observation2', 'observation3', 'observation4', 'observation5',  'observation6', 'observation7','observation8',  'observation9' )]) #which nodes were observed?
    X <- which(obs==T) #observed nodes
    Y <- round(rescaleMinMax(V(g)$reward[X], row$minList, row$maxList)) #rewards are rescaled to the same value participants saw 
    predictionNode <- row$targetNodes #which node are we trying to predict
    targetDegree <- igraph::degree(g)[predictionNode] #degree of the target node
    X.test <- seq(1:length(V(g))) #node list
    #Define kernel
    #extract parameters
    if (kernel == 'DF'){
      k <- diffusionKernel(g, alpha = alpha, normalizedLaplacian = normalizedLaplacian) 
    }else if (kernel == 'RL'){
      k <- regularizedLaplacian(g, sigma, normalizedLaplacian=normalizedLaplacian)
    }else if (kernel == 'RW'){
      k <- randomWalkKernel(g, p, normalizedLaplacian = T) #normalized laplacian needs to be the default since it isn't necessarily psd without
    }
    posterior <- gpr(X.test, X = X, Y=Y, k = k) 
    prediction[i] <-  floorCeiling(posterior[predictionNode,]$mu) # +(omegaValue * targetDegree)
    uncertainty[i] <- posterior[predictionNode,]$var
  }
  return(data.frame(mu=prediction, sig=uncertainty))
}


#dNN prediction
dNNprediction <- function(deltaValue, dataset, omegaValue=0){
  dNNpreds <- rep(NA, nrow(dataset))
  for (i in 1:nrow(dataset)){
    row <- dataset[i,]
    g <- graphList[[row$graphOrder+1]]
    maxGraphDistance =  max(distances(g)) 
    obs <- unlist(row[,c('observation1', 'observation2', 'observation3', 'observation4', 'observation5',  'observation6', 'observation7','observation8',  'observation9' )])
    X <- which(obs==T)
    Y <- NA
    Y[X] <- round(rescaleMinMax(V(g)$reward[X], row$minList, row$maxList)) #rewards are rescaled to the same value participants saw
    predictionNode <- row$targetNodes
    targetDegree <- igraph::degree(g)[predictionNode]
    distances <- distances(g)[predictionNode,] #distances between each node and prediction Node
    dnnPrediction <-  floorCeiling(mean(Y[which(distances<=deltaValue & distances>0)],na.rm=T))#average rewards of each neighbor at each distance, with floorCeiling() bounding each prediction betwteen 1-50;
    dnnPrediction[is.na(dnnPrediction)] <- 25 #if no rewards are within the distance, then predict 25
    dNNpreds[i] <- dnnPrediction + (omegaValue * targetDegree)
  }
  return(dNNpreds)
}

# kNN prediction over all k-values
kNNprediction <- function(kValue, dataset, omegaValue=0){
  kNNpreds <- rep(NA,nrow(dataset))
  for (i in 1:nrow(dataset)){
    row <- dataset[i,]
    g <- graphList[[row$graphOrder+1]]
    obs <- unlist(row[,c('observation1', 'observation2', 'observation3', 'observation4', 'observation5',  'observation6', 'observation7','observation8',  'observation9' )])
    X <- which(obs==T)
    Y <- NA
    Y[X] <- round(rescaleMinMax(V(g)$reward[X], row$minList, row$maxList)) #rewards are rescaled to the same value participants saw
    predictionNode <- row$targetNodes
    targetDegree <- igraph::degree(g)[predictionNode]
    distances <- distances(g)[predictionNode,] #distances between each node and prediction Node
    orderedNeighbors <- distances[X][order(distances[X])] #find distance of kth nearest neighbor
    orderedNeighbors <- c(orderedNeighbors, rep(max(orderedNeighbors), 7-length(orderedNeighbors))) #extend the vector to the maximum number of observed nodes using the maximum distance
    knnPrediction <- floorCeiling(mean(Y[which(distances<=orderedNeighbors[kValue] & distances>0)], na.rm=T))#average rewards of each neighbor at each distance, with floorCeiling() bounding each prediction betwteen 1-50;
    kNNpreds[i] <- knnPrediction + (omegaValue * targetDegree)
  }
  return(kNNpreds)
}


##############################################################################################################
#Standard GP Kernels
##############################################################################################################

#Radial Basis Kernel
rbf <- function(X1,X2,theta){
  #transfer to matrices
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  #check dimensions
  if(ncol(X1) != ncol(X2)){
    stop("X1 and X2 must contain input values of the same dimension.")
  } else if(!all(theta>=0)){
    stop("All parameters must be >= 0.")
  }
  #get dimensions
  N1 <- nrow(X1)
  N2 <- nrow(X2)
  d <- ncol(X1)
  #initialize sigma
  sigma <-  matrix(rep(0, N1*N2),nrow=N1)
  #observational variance
  sf <- theta[d+1]
  #noise variance
  sn <- theta[d+2]
  #loop through
  for(i in 1:d){
    #length scale
    l <- theta[i] #Note: assumes a unique length scale for each dimension
    #x-diff
    xdiff <- (outer(X1[,i],X2[,i],function(x,y) x - y)/l)^2
    sigma <- sigma + xdiff
  }
  #RBF function
  if(identical(X1,X2)){
    id <- diag(rep(1,N1))
    sigma.final <- sf*exp(-0.5*sigma) + sn*id
  } else {
    sigma.final <- sf*exp(-0.5*sigma)
  }
  #return final covariance matrix
  return(sigma.final)
}
class(rbf)<- c(class(rbf), "GP") #identify the rbf kernel as a gp model
#Ornstein-Uhlenbeck
oru <- function(X1,X2,theta){
  #matrices
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  #check dimensions
  if(ncol(X1) != ncol(X2)){
    stop("X1 and X2 must contain input values of the same dimension.")
  } else if(!all(theta>=0)){
    stop("All parameters must be >= 0.")
  }
  #get dimensions
  N1 <- nrow(X1)
  N2 <- nrow(X2)
  d <- ncol(X1) 
  #intialize matrix
  sigma <-  matrix(rep(0, N1*N2),nrow=N1) 
  #observation variance
  sf <- theta[d+1]
  #noise variance
  sn <- theta[d+2]
  #loop through
  for(i in 1:d){
    #length scale
    l <- theta[i] #Note: assumes a unique length scale for each dimension
    #x dash
    xdiff <- abs(outer(X1[,i],X2[,i],function(x,y) x - y)/l)
    sigma <- sigma + xdiff
  }
  #apply Ornstein-Uhlenbeck
  if(identical(X1,X2)){
    id <- diag(rep(1,N1))
    sigma.final <- sf*exp(-0.5*sigma) + sn*id
  } else {
    sigma.final <- sf*exp(-0.5*sigma)
  }
  #return covariance matrix
  return(sigma.final)
}

class(oru)<- c(class(oru), "GP") #identify the ornstein-uhlenbeck kernel as a gp model

#Matern 3/2
matern <- function(X1,X2,theta){
  #matrices
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  #check dimensions
  if(ncol(X1) != ncol(X2)){
    stop("X1 and X2 must contain input values of the same dimension.")
  } else if(!all(theta>=0)){
    stop("All parameters must be >= 0.")
  }
  #get dimensions
  N1 <- nrow(X1) 
  N2 <- nrow(X2) 
  d <- ncol(X1)
  #intialize matrix
  sigma <-  matrix(rep(0, N1*N2),nrow=N1)
  #observation variance
  sf <- theta[d+1]
  #noise variance
  sn <- theta[d+2]
  #loop through
  for(i in 1:d){
    #length scale
    l <- theta[i] #Note: assumes a unique length scale for each dimension
    #x dash
    xdiff <- abs(outer(X1[,i],X2[,i],function(x,y) x - y)/l)
    sigma <- sigma + xdiff
  }
  #apply Mater 3/2
  if(identical(X1,X2)){
    id <- diag(rep(1,N1))
    sigma.final <- sf*(1+(sqrt(3)*sigma))*exp(-sqrt(3)*sigma) + sn*id
  } else {
    sigma.final <- sf*(1+(sqrt(3)*sigma))*exp(-sqrt(3)*sigma)
  }
  #return covariance matrix
  return(sigma.final)
}

class(matern)<- c(class(matern), "GP") #identify the rbf kernel as a gp model

##############################################################################################################
#MATRIX INVERSION
##############################################################################################################

#calculate inverse of the cov-function using sigular value decomposition
cov.inverse.svd <- function(X, tol = sqrt(.Machine$double.eps)){
  # Generalized Inverse of a Matrix
  dnx <- dimnames(X)
  if(is.null(dnx)) dnx <- vector("list", 2)
  #singular value decomposition
  s <- svd(X)
  nz <- s$d > tol * s$d[1]
  #inverse
  K.inv <- structure(
    if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X,
    dimnames = dnx[2:1])
  #logarithm of determinant.
  #log.K.det <- sum(log(s$d))
  return(K.inv)
}

#calculate inverse of the cov-function using Cholesky
cov.inverse.chol <- function(X){
  #cholseky decomposition
  R <- chol(X)
  #complex conjugate
  Rt <- Conj(t(R))
  #invert
  R.inv <- solve(R)
  #invert
  Rt.inv <- solve(Rt)
  #multiply matrices
  X.inv <- R.inv %*% Rt.inv
  #log determinant
  #log.X.det <- 2*sum(log(diag(R))) 
  return(X.inv)
}

##############################################################################################################
#GAUSSIAN PROCESS
##############################################################################################################

#Gaussian Process function
#X.test: matrix for predcitions
#X; matrix of observations
#y: vector of observed outcomes
#k: covariance matrix (e.g., from diffusion kernel)
gpr<- function(X.test, X, Y, k, mu_0 = 25, noise = 0.0001){
  Y <- Y - mu_0 #invariant mean scaling
  #make it a matrix
  Xstar <- as.matrix(X.test)
  #dimensions
  d <- ncol(as.matrix(X))
  #Calculate K:
  K <- as.matrix(k[X,X]) +  (diag(length(Y))*noise) 
  KK.inv <- chol2inv(chol(K)) #MASS implementation of Cholesky  
  #times y
  Ky <- KK.inv %*% Y
  #apply the kernel
  result <- apply(Xstar, 1, function(x){
    XX <- matrix(x,nrow=1) 
    Kstar <- k[X, XX]
    Kstarstar <- k[XX,XX]
    #get mean vector
    mu <- t(Kstar) %*% Ky
    #get covariance
    cv <- Kstarstar - (t(Kstar) %*% KK.inv %*% Kstar) 
    if (cv<0){ cv <- abs(cv)} #DEBUG SOLUTION: MANUALLY SET CV TO BE POSITIVE IF NEGATIVE
    #return means and variance
    return(c(mu, cv))
  })
  #as a data frame with names mu and sig
  prediction <- as.data.frame(t(result))
  colnames(prediction) <- c("mu", "var")
  prediction$mu <- prediction$mu + mu_0 #account for mean function
  #image(matrix(prediction$mu, nrow=10, ncol = 10)) #mean reward
  #image(matrix(prediction$sig, nrow=10, ncol = 10)) #mean uncertainty
  return(prediction)
}

##############################################################################################################
## KALMAN FILTER
##############################################################################################################

#Like a kalman filter, but with stable-mean and no innovation variance
bayesianMeanTracker <- function(x, y, kNoise=1, prevPost=NULL, uniqueOptions = 64, mu0=0){ #Updates the previous posterior based on a single observation
  #parameters
  var0 <- 5 #prior variance
  vare <- kNoise #error varriance
  if (is.null(prevPost)){#if no posterior prior, assume it is the first observation
    predictions <- data.frame(mu=rep(mu0,uniqueOptions), var=rep(var0,uniqueOptions))
  }else{#if previous posterior is provided, update
    predictions <- prevPost
  }
  #Which of the 64 options were chosen at time?
  allopts<- 1:uniqueOptions
  chosen <- x
  #Kalman gain
  kGain <- predictions$var[chosen] / (predictions$var[chosen] + vare^2)
  #update mean
  predictions$mu[chosen] <- predictions$mu[chosen] + (kGain * (y-predictions$mu[chosen]))
  #update variance for observed arm
  predictions$var[chosen] <- predictions$var[chosen] * (1 - kGain)
  #return output
  return(predictions)
}
class(bayesianMeanTracker)<- c(class(bayesianMeanTracker), "KalmanFilter")


##############################################################################################################
#SIMPLE PREDICTION MODELS
##############################################################################################################
gpPred <- function(g, target, betaValue){
  k <- diffusionKernel(g, beta = betaValue)
  X <- which(!is.na(V(g)$reward))
  Y <- as.numeric(V(g)[X]$reward)
  return(gpr(V(g), X, Y, k)[target,"mu"])
}

nnPred <- function(g, target, dist){
  distances <- distances(g)[target,] #distances between each node and prediction Node
  X <- which(!is.na(V(g)$reward))
  Y <- as.numeric(V(g)[X]$reward)
  dnnPrediction <- floorCeiling(mean(Y[which(X %in% which(distances<=dist & distances>0))]))#average rewards of each neighbor at each distance, with floorCeiling() bounding each prediction betwteen 1-50; 
  dnnPrediction[is.na(dnnPrediction)] <-25 #if no rewards are within the distance, then predict 25
  return(dnnPrediction)
}



##############################################################################################################
#SIMPLE Choice MODELS
##############################################################################################################

dnnChoice <- function(distanceMatrix, X, Y, dist, mu_0 = 0 ){
  values <- rep(NA, 64)
  values[X] <- Y
  dnnPrediction <- sapply(1:64, FUN=function(j)  mean(values[which(distanceMatrix[j,]<=dist & distanceMatrix[j,] >0)],na.rm=T)) #average rewards of each neighbor within distance
  dnnPrediction[is.na(dnnPrediction)] <-mu_0 #if no rewards are within the distance, then predict mu_0
  dnnPrediction[X] <- Y #default predictions of previously observed nodes to the last reward observed
  prediction <- data.frame(mu = dnnPrediction, var = rep(0,64))
  return(prediction)
}

class(dnnChoice)<- c(class(dnnChoice), "dNN")

knnChoice <- function(distanceMatrix, X, Y, kValue, mu_0 = 0 ){
  values <- rep(NA, 64)
  values[X] <- Y
  knnPrediction <- sapply(1:64, FUN=function(j) {
    orderedNeighbors <- distanceMatrix[j,][order(distanceMatrix[j,])] #find distance of kth nearest neighbor
    maxDistance <- orderedNeighbors[kValue+1]
    mean(values[which(distanceMatrix[j,]<=maxDistance & distanceMatrix[j,] >0)], na.rm=T)
  }) #average rewards of each neighbor within distance
  knnPrediction[is.na(knnPrediction)] <-mu_0 #if no rewards are within the distance, then predict mu_0
  knnPrediction[X] <- Y #default predictions of previously observed nodes to the last reward observed
  prediction <- data.frame(mu = knnPrediction, var = rep(0,64))
  return(prediction)
}

class(knnChoice)<- c(class(knnChoice), "kNN")

nullChoice <- function(){ #Null model used for stickiness only
  return (data.frame(mu=rep(0,64), var = rep(0,64)))
}

class(nullChoice) <- c(class(nullChoice), 'null')

##############################################################################################################
#ACQUISITION FUNCTIONS
##############################################################################################################

ucb<- function(means, variances, beta){
  return (means + (beta *sqrt(variances) ))
}
class(ucb)<- c(class(ucb), "UCB")

#Probability of Maximum Utility (the model formerly known as Thompson Sampling)
pmu<-function(out){
  #initialize
  mus<-matrix(out$mu, nrow=nrow(out)/uniqueOptions, byrow=TRUE)
  sigs<-matrix(sqrt(out$sig), nrow=nrow(out)/uniqueOptions, byrow=TRUE) #convert variance to standard deviation
  outtotal<-mus
  #loop through
  for (i in 1:nrow(mus)){
    #simulate from normals given the means and standard deviation
    simout<-apply(cbind(mus[i,], sigs[i,]),1, function(x){rnorm(n=1000,mean=x[1],sd=x[2])})
    #collect mean simulated wins per arm
    outtotal[i,]<-apply(simout==apply(simout,1,max),2,mean) 
  }
  #return them
  return(outtotal+0.0001)
}

#Probability of Improvement
poi<-function(out, y.star){
  #calulate the probabilities of improvement
  #pi<-pnorm((out$mu-y.star)/out$sig)
  laplace <- sqrt(.Machine$double.eps) #DEBUG SOLUTION TO NAN VALUES; LAPLACIAN SMOOTHING
  pi<-pnorm((out$mu-y.star+laplace)/(sqrt(out$sig)+laplace)) 
  outtotal<-matrix(pi, nrow(out)/uniqueOptions, byrow=TRUE)
  #return them
  return(outtotal)
}
#add "IMP" to the class of the probofimp function, so that modelFit can recognize that it takes y.star as an argument
class(poi)<- c(class(poi), "Imp")

#Expected improvement
exofimp<-function(out, y.star){
  #calulate z-scores first, then expected improvements
  laplace <- sqrt(.Machine$double.eps) #DEBUG SOLUTION TO NAN VALUES; LAPLACIAN SMOOTHING
  z<-(out$mu-y.star+laplace)/(sqrt(out$sig)+laplace)
  ei <-(out$mu-y.star)*pnorm(z)+sqrt(out$sig)*dnorm(z)
  outtotal<-matrix(ei, nrow(out)/uniqueOptions, byrow=TRUE)
  #return them
  return(outtotal)
}
#add "IMP" to the class of the exofimp function, so that modelFit can recognize that it takes y.star as an argument
class(exofimp)<- c(class(exofimp), "Imp")

#Greedy Mean
greedyMean <- function(out, uniqueOptions=64){
  outtotal<-out$mu #the value of each choice is solely based on the expectation of reward
  #avoid borderline cases
  #outtotal[outtotal<0]<-0.0001
  #outtotal[outtotal>100]<-100
  outtotal<-matrix(outtotal, nrow=nrow(out)/uniqueOptions, byrow=TRUE)
}
class(greedyMean)<- c(class(greedyMean), "greedy")

#Greedy Variance
greedyVar <- function(out, uniqueOptions=64){
  outtotal<- sqrt(out$var) #the value of each choice is solely based on the expected uncertainty
  #avoid borderline cases
  #outtotal[outtotal<0]<-0.0001
  #outtotal[outtotal>100]<-100
  outtotal<-matrix(outtotal, nrow=nrow(out)/uniqueOptions, byrow=TRUE)
}
class(greedyVar)<- c(class(greedyVar), "greedy")


########################################################################################################################
# Parameter labeller
########################################################################################################################

paramLabeller <- function(k, samplingStrategy, stickiness=F){
  params <- c()
  #Sampling strategies first
  if (inherits(samplingStrategy, "UCB")){
    params <- c(params, c('tau', 'beta'))
  }else if(inherits(samplingStrategy, "greedy")){
    params <- c(params,'tau')
  }
  #learning models
  if (inherits(k, "RBF") ){
    params <- c(params, 'lambda')  
  }else if (inherits(k, "DF")){
    params <- c(params, 'alpha')  
  }else if (inherits(k, "KalmanFilter")){
    params <- c(params, 'errorVariance')  
  }else if (inherits(k, "dNN")){
    params <- c(params, 'd')
  }else if (inherits(k, 'kNN')){
    params <- c(params, 'k')
  }
  #stickiness
  if (stickiness){
    params <- c(params, 'omega')
  }
  return(params)
}
