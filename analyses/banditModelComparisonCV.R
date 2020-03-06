#Model comparison
rm(list=ls()) #house keeping

#load packages
packages <- c('plyr', 'jsonlite', 'DEoptim', 'igraph', "matrixcalc", "fields", "MASS")
lapply(packages, require, character.only = TRUE)

source('exportImportGraph.R')
source('models.R')
source('utilities.R')
##############################################################################################################
#Cluster configuration: (1 subject x model combination) per CPU
##############################################################################################################

#IMPORTANT: update batch name
batchName = 'networkBandit' #saves output csv files in 'modelResults/batchName/*.csv'
sampleSize <- 98
stickiness = T #If this is turned off, then also change the parameter mapping function for dNN and kNN 

#create list of all kernel functions (i..e, learning rules)
kernellist<-list(diffusionKernel, bayesianMeanTracker, knnChoice, dnnChoice, nullChoice) #GP, BMT, kNN, dNN, and null choices (i.e., stickiness only model)


#names of all kernel functions
kernelnames<-c('DF', 'BMT', 'kNN', 'dNN',  'Sticky')


#list of all acquisition functions
acqlist<-list(greedyVar, greedyMean, ucb) #different sampling strategies. Only UCB reported, although greedy variance and greedy mean are lesioned versions of UCB

#names of all acquisition functions
acqnames<-c("GV", "GM", "UCB")

#Cluster id from qsub
#clusterid <- sample(1:(length(kernellist)*length(acqlist)*sampleSize),1)
clusterid <- as.integer(commandArgs(TRUE)[1]) #Cluster id, corresponds to an integer used to indicate which combination of kernel and acquisition function to simulate

#all combinations of kernels and acquisition functions will be needed
combs<-expand.grid(1:length(kernellist), 1:length(acqlist))

#create a matrix with  combinations of subjectIds and model combinations
subjectComb <- expand.grid(1:sampleSize, 1:(length(kernellist) * length(acqlist))) #1:? defines the number of unique models to be analyzed
subjectId <- subjectComb[clusterid,1] #used to identify unique subjects
combId <- subjectComb[clusterid,2] #used to identify a unique model (i.e., combination of GP kernel and acquisition function)

#trim combs if cluster==TRUE
model <- combs[combId,]

set.seed(clusterid) #set seed as the clusterid

#Parameters
#parameter bounds; all parameters  (except d and k) are exponentiated inside the optimization function
paramUBound <- list(tau = 5, beta = 5, lambda = 5, errorVariance  = 5, alpha = 2, d = 20, k = 25,  omega = 5) 
paramLBound <- list(tau = -5,  beta = -5, lambda = -5, errorVariance  = -5, alpha = -4, d = 0, k =0, omega = -5 ) 

params <- paramLabeller(kernellist[[model[[1]]]], acqlist[[model[[2]]]],  stickiness)
allParams <- rep(NA,9)
names(allParams) <- c('tau', 'beta', 'lambda', 'errorVariance', 'alpha', 'd', 'k', 'omega') #all parameters
##############################################################################################################
#Compile Experimental Data
##############################################################################################################
data <- dataImportBandit('../data/banditTask/networkBandit.csv', delim=',')
data <- subset(data, round<10) #remove bonus round
#Normalize data
data$z <- (data$rewardObjective - 50) / 100 #convert to zero mean and unit variance
uid <- unique(data$id)[subjectId] #convert subjectId to uid


####################################################################################
# Load pre-generated graphs
####################################################################################
experimentGraphFiles <- '../data/banditTask/network2.json'
networks <- fromJSON(experimentGraphFiles)
edges <- networks[[1]]
nodes <- networks[[2]]

graphList = list()
distanceTables = list()
nodeDegrees = list()

for (i in 1:length(edges)){
  g <- graph_from_edgelist(as.matrix(edges[[i]]), directed=F) 
  V(g)$reward  <- nodes[[i]]$reward
  #plot(g, vertex.size = V(g)$reward)
  graphList[[i]] <-g
  distanceTables[[i]] <- distances(g)
  nodeDegrees[[i]] <- igraph::degree(g)
}

##############################################################################################################
#Model Fitting
##############################################################################################################

#Negative Log Likelihood of model prediction 
#parameters are distinct for each individual model
#subjD = subset of data for a specific subject
#kernel = any of fastKalmanFilter, standard GP (rbf, matern, or oru), or lingp
#acquisition = any of ucb, probofimp, expofimp, or pmu
#rounds list of numbers between 1:4, relative to the order of rounds with the specified horizon length
#refactorUCB: whether or not to refactor tau and beta parameters in UCB
#returnPredictions: whether or not to also return the dataframe of probabilistic model predictions for each choice. False for MLE, True for out of sample prediction
modelFit<-function(par, subjD, acquisition, k, rounds, returnPredictions=F){
  #Extract and process parameters
  names(par) <- params
  tau<-exp(par['tau'])  #inverse temperature for softmax; all models have it
  #Which posterior function to use; therefore, which parameters to use
  if (inherits(k, "KalmanFilter")){ #kalman filter model
    kNoise <- exp(par['errorVariance'])
  }else if(inherits(k, "RBF")){ #RBF kernel uses lambda
    lambda <- exp(par['lambda'])
    K <- rbf(gridLoc, gridLoc, lambda = lambda) #compute gram Matrix
  }else if(inherits(k, "DF")){ #Diffusion kernel uses alpha
    alpha <- exp(par['alpha'])
  }else if(inherits(k, 'dNN')){
    dValue <- round(par['d']) #make sure it is an integer
  }else if(inherits(k, 'kNN')){
    kValue <- round(par['k']) #make sure it is an integer
  }
  #Add stickiness parameter 
  if (stickiness==T){
    omega <- exp(par['omega'])
  }
  #Additional acquisition function dependent parameters
  if (inherits(acquisition, "UCB")){ #check if UCB is used
    beta <- exp(par['beta']) #If UCB, beta is always 2nd last
  }
  #Vector to store negative log likelihods
  nLL <- rep(0,length(rounds)) 
  for (r in rounds){ #Begin looping through each round
    #subset of data for round r
    roundD <- subset(subjD, round==r)
    #Observations of subject choice behavior
    Observed <- roundD[1:25, 'nodeId']
    chosen <- roundD[2:26,'nodeId'] 
    y  <- roundD$z[1:25] #trim off the last observation, because it was not used to inform a choice (round already over)
    g<- graphList[[roundD$graphOrder[[1]]]] #extract graph
    distanceMatrix <- distanceTables[[roundD$graphOrder[[1]]]] 
    #Compute diffusion kernel gram matrix
    if(inherits(k, "DF")){ #Diffusion kernel 
      K <- diffusionKernel(g,alpha=alpha)
    }
    #Utilties of each choice
    utilities <-NULL
    prevPost <- NULL
    #loop through observations
    for (i in 1:25){ #skip the first observation, since that was completely random
      #Which posterior function to use
      if (inherits(k, "KalmanFilter")){# kalman filter/Bayesian mean tracker model
        out<- bayesianMeanTracker(x = Observed[i], y=y[i], prevPost = prevPost, kNoise = kNoise)
        prevPost <- out #update prevPost for the next round
      }else if (inherits(k, 'dNN')){ #d nearest neighbors models
        out <- dnnChoice(distanceMatrix, X=Observed[1:i], Y=y[1:i], dist = dValue, mu_0=0)
      }else if (inherits(k, 'kNN')){ #k nearest neighbors models
        out <- knnChoice(distanceMatrix, X=Observed[1:i], Y=y[1:i],kValue=kValue, mu_0=0)
      }else if (inherits(k,'null')){
        out <- nullChoice()
      }else{ #Both RBF and Diffusion kernel implementations of the GP use the same function, based on the gram matrix K
        out <- gpr(X.test=1:64, X=Observed[1:i], Y=y[1:i], k=K, mu_0=0)
      }
      #Slightly different function calls for each acquisition function
      if (inherits(acquisition, "UCB")){ #UCB takes a beta parameter
        utilityVec<-ucb(out$mu, out$var, beta)
      }else if(inherits(acquisition, "Imp")){ #ProbImp or ExpectedImp
        y.star <- max(roundD$z[1:i]) #Best revealed solution
        utilityVec <- acquisition(out, y.star)
      }else{ #PMU or any other greedy model
        utilityVec <- acquisition(out)
      }
      #Add stickiness if specified
      if (stickiness==TRUE){ 
        #weight the utilityVec by the inverse manhattan distance
        utilityVec[Observed[i]] <- utilityVec[Observed[i]] + omega
      }
      utilityVec <- utilityVec - max(utilityVec) #avoid overflow
      utilities <- rbind(utilities, utilityVec) # build horizon_length x 25 matrix, where each row holds the utilities of each choice at each decision time in the search horizon
    }
    #Softmax rule
    p <- exp(utilities/tau)
    p <- p/rowSums(p)
    #avoid underflow by setting a floor and a ceiling
    p <- (pmax(p, 0.00001))
    p <- (pmin(p, 0.99999))
    #Calculate Negative log likelihood
    nLL[which(r==rounds)] <- -sum(log(p[cbind(c(1:25),chosen)]))
  }#end loop through rounds
  if (returnPredictions==FALSE){ #Return only summed log loss for computing MLE
    return(sum(nLL))  #Return negative log likelihoods of all observations 
  }else if (returnPredictions==TRUE){ #Return detailed dataframe of model predictions for outofsample predictions
    detailedPredictions <- list(sum(nLL), p, roundD[2:26,'nodeId']) #construct a list of summed log loss, model predictions, and chosen item
    return(detailedPredictions) #return detailed predictions
  }
}

##############################################################################################################
#CROSS VALIDATION FUNCTION
##############################################################################################################
#function to plug in to the optimaztion routine
#selector: scalar, indicates a specific participant
#kernel function
#acquisition, function,
#leaveoutindex is the round being held out
cvfun<-function(selector, kernelFun, acquisition, leaveoutindex){
  #subselect participant, horizon and rounds not left out
  d1<-subset(data, id==selector)
  #training set
  rounds <- seq(1,9) #Trim last round, since bonus round isn't modeled
  trainingSet <- rounds[! rounds==rounds[leaveoutindex]] #remove round specified by leaveoutindex
  #test set
  testSet <- rounds[leaveoutindex]
  #Compute MLE on the trainingSet
  lbound<-unlist(paramLBound[params])
  ubound<-unlist(paramUBound[params])
  #Begin cross validation routine
  if (length(params)>=2){#if 2 or more parameters
    #TRAINING SET
    #enforce cardinality constraints
    if (inherits(kernelFun,'dNN') | inherits(kernelFun, 'kNN')){
      fnmap_f <- function(x) c(x[1],round(x[2]),x[3]) #round the k and d parameters, since they are an integer
      fit<-DEoptim(modelFit, lower=lbound, upper=ubound, subjD=d1, k=kernelFun, rounds = trainingSet, acquisition=acquisition, DEoptim.control(itermax=100), fnMap = fnmap_f)
    }else{ #normal case
      fit<-DEoptim(modelFit, lower=lbound, upper=ubound, subjD=d1, k=kernelFun, rounds = trainingSet, acquisition=acquisition,DEoptim.control(itermax=100))
    }
    paramEstimates <- fit$optim$bestmem #MODEL DEPENDENT PARAMETER ESTIMATES
    allParams[params] <- paramEstimates #index all possible parameters by the parameters of this given model, and insert those values
    #TEST SET
    predict <- modelFit(par=paramEstimates, subjD=d1, acquisition=acquisition, k=kernelFun, rounds=testSet,  returnPredictions=TRUE)
    cvresults <- data.frame(loo = leaveoutindex,nLL =  predict[[1]], tau = exp(allParams['tau']), beta = exp(allParams['beta']), lambda = exp(allParams['lambda']), errorVariance = exp(allParams['errorVariance']), alpha = exp(allParams['alpha']), d = round(allParams['d']), k = round(allParams['k']),  omega = exp(allParams['omega'])) #leaveoutindex, nLL, parameters....
    output <- list(cvresults, predict[[2]], predict[[3]])
  } else{
    #TRAINING SET
    fit<-optimize(modelFit, lower=lbound, upper=ubound, subjD=d1, k=kernelFun, rounds = trainingSet, acquisition=acquisition)
    paramEstimates <- fit$minimum #MODEL DEPENDENT PARAMETER ESTIMATES
    allParams[params] <- paramEstimates #index all possible parameters by the parameters of this given model, and insert those values
    #TEST SET
    predict <- modelFit(par=paramEstimates, subjD=d1, acquisition=acquisition, k=kernelFun, rounds=testSet, returnPredictions=TRUE)
    cvresults <- data.frame(loo = leaveoutindex,nLL =  predict[[1]], tau = exp(allParams['tau']), beta = exp(allParams['beta']), lambda = exp(allParams['lambda']), errorVariance = exp(allParams['errorVariance']),  alpha = exp(allParams['alpha'])) #leaveoutindex, nLL, parameters....
    output <- list(cvresults, predict[[2]], predict[[3]])
  }
  return(output) #return optimized value
}


##############################################################################################################
#OPTIMIZATION ROUTINE
##############################################################################################################


#cross-validation routine
crossvalidation <- data.frame() #cross-validation results
modelPredictions <- lapply(1:9, matrix, data= NA, nrow=25, ncol=64)
chosenMatrix <- matrix(NA, nrow=9, ncol=25)


start.time <- Sys.time()
for (loo in 1:9){#leave one out index
  cv <- cvfun(selector = uid, kernelFun=kernellist[[model[[1]]]], acquisition = acqlist[[model[[2]]]], leaveoutindex=loo)
  crossvalidation <- rbind(crossvalidation, cv[[1]])
  modelPredictions[[loo]] <- cv[[2]]
  chosenMatrix[loo,] <- cv[[3]]
}

output <- list(crossvalidation, modelPredictions, chosenMatrix)
#save the vector with kernel-acquisition-pair as name
name<-paste0("modelResults/", batchName, "/", kernelnames[model[[1]]], acqnames[model[[2]]], uid)
save(output, file=paste0(name, '.Rdata'))

print(output)

end.time <- Sys.time()
elapsed <- end.time - start.time
print(elapsed)

