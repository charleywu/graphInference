#Exp1 modeling
rm(list=ls())

packages <- c( 'tidyr',"jsonlite", "igraph", 'matrixcalc', 'Matrix', 'dplyr', 'DEoptim')
lapply(packages, require, character.only = TRUE)

source('exportImportGraph.R')
source('models.R')
source('utilities.R')
#source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

batchname <- 'fullPaper'

# clusterid <- round(runif(1,1,100))
clusterid <- as.integer(commandArgs(TRUE)[1]) #for MPI cluster
set.seed(clusterid)
####################################################################################
# Load participant data and graphs
####################################################################################

df <- dataImportExp1('../data/exp1/experiment1Full.csv', delim=',')

experimentGraphFiles <- '../data/exp1/graphs.json'
networks <- fromJSON(experimentGraphFiles)
edges <- networks[[1]]
nodes <- networks[[2]]

graphList = list()

for (i in 1:length(edges)){
  g <- graph_from_edgelist(as.matrix(edges[[i]]), directed=F) 
  V(g)$reward  <- nodes[[i]]$reward
  #plot(g, vertex.size = V(g)$reward)
  graphList[[i]] <-g
}



####################################################################################
# DEFINE LOSS FUNCTIONS
#
####################################################################################

#Loss Function for Gaussian Process predictions 
gpErrorFunc <- function(pars, kernel, dataset, normalize = F){ #This mainly just parses two parameters if necessary
  if (kernel == 'DF'){
    alpha = pars[1]
  }else if (kernel == 'RL'){
    sigma = pars[1]
  }else if (kernel == 'RW'){
    p = pars[1]
  }
  modelPreds <-  gpPrediction(pars,kernel=kernel, dataset = dataset, normalizedLaplacian=normalize)$mu
  humanPreds <- dataset$predictions
  SE <- (modelPreds - humanPreds)^2 #squared error
  SSE <- sum(SE) #sum of squared errors
  return(SSE) #return sum of squared errors
}


#d-Nearest Neighbors with omega
dNNErrorFunc <- function(deltaValue,  dataset){ #Delta is the maximum distance at which NN includes nodes to be averaged
  errorVec <-rep(NA, nrow(dataset))
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
    dnnPrediction <- floorCeiling(mean(Y[which(distances<=deltaValue & distances>0)],na.rm=T))#average rewards of each neighbor at each distance, with floorCeiling() bounding each prediction betwteen 1-50;
    dnnPrediction[is.na(dnnPrediction)] <- 25 #if no rewards are within the distance, then predict 25
    dnnPrediction <- dnnPrediction #+ (omegaValue*targetDegree)
    #squared error
    errorVec[i] <- (dnnPrediction - row$trueTargetValue)^2
  }
  return(sum(errorVec)) #return sum of squares  
}


#k-Nearest Neighbors with omega
kNNErrorFunc <- function(kValue,  dataset){ #Delta is the maximum distance at which NN includes nodes to be averaged
  errorVec <-rep(NA, nrow(dataset))
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
    orderedNeighbors <- distances[X][order(distances[X])] #find distance of kth nearest neighbor
    orderedNeighbors <- c(orderedNeighbors, rep(max(orderedNeighbors), 7-length(orderedNeighbors))) #extend the vector to the maximum number of observed nodes using the maximum distance
    knnPrediction <- floorCeiling(mean(Y[which(distances<=orderedNeighbors[kValue] & distances>0)], na.rm=T))#average rewards of each neighbor at each distance, with floorCeiling() bounding each prediction betwteen 1-50;
    knnPrediction <- knnPrediction  #+ (omegaValue*targetDegree)
    #squared error
    errorVec[i] <- (knnPrediction - row$trueTargetValue)^2
  }
  return(sum(errorVec)) #return sum of squares  
}



####################################################################################
# Model fitting
#
####################################################################################

k_id <- unique(df$id)[clusterid] #each iteration is one participant
subdf <- subset(df, id==k_id)

#add to dataframe
newColumns <- c('alpha', 'gpDFpred', 'gpDFuncertainty', 'delta','dNNpred', 'k','kNNpred')
df[,newColumns] <- NA
#Define integer constraints on Random Walk kernel p parameter
fnmap_RW <- function(x)c(round(x))

start <- proc.time()

#1. Simple models first
#standard dNN
dNNpreds <- sapply(1:6, function(d) dNNprediction(d,  dataset=subdf)) #Make predictions
dnnError <- (dNNpreds-subdf$predictions)^2 #compute squared error vs. participant predictions
dNN_mle <- t(sapply(1:nrow(subdf), function(i) {
  test <- dnnError[i,]
  training <- dnnError[-i,]
  delta <- sample(which(colSums(training)==min(colSums(training))),1) #choose min, and break ties at random
  return(c(delta,dNNpreds[i,delta]))
}))

#Put in the values
df[as.numeric(rownames(subdf)),'delta'] <-dNN_mle[,1] #which value had the lowest MSE is the each training set
df[as.numeric(rownames(subdf)),'dNNpred'] <- dNN_mle[,2]

#standard kNN
kNNpreds <- sapply(1:7, function(k) kNNprediction(k, dataset= subdf))
knnError <- (kNNpreds-subdf$predictions)^2
kNN_mle <- t(sapply(1:nrow(subdf), function(i) {
  test <- knnError[i,]
  training <- knnError[-i,]
  kValue <- sample(which(colSums(training)==min(colSums(training))),1) #choose min, and break ties at random
  c(kValue,kNNpreds[i,kValue])
}))
#Put in the values
df[as.numeric(rownames(subdf)),'k'] <-kNN_mle[,1]
df[as.numeric(rownames(subdf)),'kNNpred'] <-kNN_mle[,2]

#2. Fit GP which require optimization
for (i in 1:nrow(subdf)){ #loop through rounds
  testSet <- subdf[i,]
  trainingSet <- subdf[-i,]
  #Fit GP-Diffusion kernel
  gpDF_mle <- optimize(f = gpErrorFunc, interval =c(0.0001, 4), kernel = 'DF', dataset = trainingSet)
  df[as.numeric(rownames(subdf)[i]), 'alpha'] <- gpDF_mle$minimum
  preds <- gpPrediction(gpDF_mle$minimum, kernel = 'DF', dataset=testSet)
  df[as.numeric(rownames(subdf)[i]), 'gpDFpred'] <- preds$mu
  df[as.numeric(rownames(subdf)[i]), 'gpDFuncertainty'] <- preds$sig
}

end <-  proc.time()
print(end-start)

#Save data
saveRDS(df[as.numeric(rownames(subdf)),] , file =paste0('modelResults/Exp1/',batchname,'/',clusterid,'.RDS'))

