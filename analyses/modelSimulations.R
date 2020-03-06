#Charley Wu 2019

#Script to run rational models in comparison to human behavior

#############################################################################################################################
# IMPORT DATA AND ENVIRONMENTS
#############################################################################################################################
#house keeping
rm(list=ls())

#load packages
packages <- c('tidyr',"jsonlite", 'plyr', "igraph", 'matrixcalc', 'Matrix', 'dplyr', 'stargazer', 'coefplot', "grid", 'matrixcalc', 'parallel')
lapply(packages, require, character.only = TRUE)

source('exportImportGraph.R')
source('models.R')
source('utilities.R')

reps <- 10000 #replications
cores <- detectCores()-2

##############################################################################################################
#Load pre-generated graphs
##############################################################################################################
experimentGraphFiles <-'../data/banditTask/network2.json'
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

#############################################################################################################################
# RANDOM MODEL
# to read the previous simulation from disk:
# randomDF <- read.csv("rationalModels/random.csv")
#############################################################################################################################

randomModel <- function(replications, outputfile){
  #Run for both types of environments
  rewards <- mcmapply(1:replications, FUN=function(x) sample(V(graphList[[sample(1:40,1)]])$reward, 26, replace = TRUE), mc.preschedule=TRUE, mc.cores=cores)
  #put into dataFrame
  randomDF <- data.frame(trial=seq(1:26)-1,
                         meanReward = rowMeans(rewards), #reward averaged over replications
                         meanSE = apply(rewards, 1, FUN = function(x) sd(x)/sqrt(length(x))))
  #add model label
  randomDF$Model <- "Random"
  #write output
  if(!is.null(outputfile)){
    write.csv(randomDF, outputfile) 
  }
  return(randomDF)
}

randomDF <- randomModel(reps, "rationalModels/random.csv")

#############################################################################################################################
# BMT-UCB Model
# to read the previous simulation from disk:
# bmtDF <- read.csv("rationalModels/BMTUCB.csv")
#############################################################################################################################

bmtRationalModel <- function(replications, modelName, outputfile, parameters, stickiness=T){
  #run for smooth environments
  rewardTotal <- mcmapply(1:replications, FUN=function(x){
    #sample parameters
    params <- sample_n(parameters,1)
    kError <- params$errorVariance
    beta <- params$beta
    tau <- params$tau
    omega <- params$omega
    #tau <- median(smoothPars$tau)
    #randomly choose environment
    envNum <- sample(1:40,1) 
    #1st trial is random
    location <- sample(1:64,1)#each location as an integer from 1:25; first value is random, so treat both environments as replications of the same actions
    #Observations
    g <- graphList[[envNum]]
    #reward
    reward<- c()
    reward[1] <- Y <- V(g)$reward[location] #first sample is random
    #posterior
    prevPost <- NULL  #initialize previous posterior for sequential updating of BMT posterior
    for (j in 2:26){ #after that, loop through remaining trials and make decisions based on BMT-UCB preditions
      #update posterior predictions
      post <- bayesianMeanTracker(x =location, y = (V(g)$reward[location] - 50)/100, prevPost = prevPost, kNoise = kError) #normalize observation of reward
      prevPost <- post  #save new posterior as prevPost for next round
      #compute acquisition function evaluation
      utilityVec <- ucb(post$mu,post$var, beta = beta)
      if (stickiness==T){ #add stickiness
        utilityVec[location]<- utilityVec[location] + omega
      }
      #to prevent overflow, subtract max of q(x) vector 
      utilityVec <- utilityVec - max(utilityVec)
      #compute softmax choice probabilities
      p <- exp(utilityVec/tau)
      p <- p/sum(p)
      #Sample next choice
      location <- sample(1:64,1, prob=p, replace=TRUE)
      #update reward, X1, X2, and Y for both smooth and rough
      reward[j] <- V(g)$reward[location]
      }
    reward
    }, mc.preschedule = TRUE, mc.cores=cores)
  
  #put into dataFrame
  bmtDF <- data.frame(trial=seq(1:26)-1,
                     meanReward = rowMeans(rewardTotal), #reward averaged over replications
                     meanSE = apply(rewardTotal, 1, FUN = function(x) sd(x)/sqrt(length(x))))
  #add model label
  bmtDF$Model <- modelName
  #write to csv
  if(!is.null(outputfile)){
    write.csv(bmtDF, outputfile)
  }
  return(bmtDF)
}

bmtPars <- read.csv('rationalModels/parameters/BMT.csv')
bmtDF <- bmtRationalModel(reps, modelName = 'BMT-UCB', outputfile = "rationalModels/BMTUCB.csv", parameters = bmtPars)
#############################################################################################################################
# GP Model
# to read the previous simulation from disk:
# gpDF <- read.csv("rationalModels/GPUCB.csv")
#############################################################################################################################

gpRationalModel <- function(replications, modelName, outputfile, parameters, stickiness = T){
  #run for smooth environments
  rewardTotal <- mcmapply(1:replications, FUN=function(x){
    #sample parameters
    params <- sample_n(parameters,1)
    alpha <- params$alpha
    beta <- params$beta
    tau <- params$tau
    omega <- params$omega
    #tau <- median(smoothPars$tau)
    #randomly choose environment
    envNum <- sample(1:40,1) 
    observations <- rep(NA,26)
    #1st trial is random
    location <-observations[1]<- sample(1:64,1)#each location as an integer from 1:25; first value is random, so treat both environments as replications of the same actions
    #Observations
    g <- graphList[[envNum]]
    K <- diffusionKernel(g, alpha)
    #reward
    reward<- rep(NA, 26)
    reward[1] <- V(g)$reward[location] #first sample is random
    for (j in 2:26){ #after that, loop through remaining trials and make decisions based on BMT-UCB preditions
      #update posterior predictions
      post <- gpr(seq(1:64),X = observations[1:(j-1)], Y = (reward[1:(j-1)]- 50)/100, k=K, mu_0=0) #normalize observation of reward
      #compute acquisition function evaluation
      utilityVec <- ucb(post$mu,post$var, beta = beta)
      if (stickiness==T){ #add stickiness
        utilityVec[location]<- utilityVec[location] + omega
      }
      #to prevent overflow, subtract max of q(x) vector 
      utilityVec <- utilityVec - max(utilityVec)
      #compute softmax choice probabilities
      p <- exp(utilityVec/tau)
      p <- p/sum(p)
      #Sample next choice
      observations[j] <-  sample(1:64,1, prob=p, replace=TRUE)
      #update reward, X1, X2, and Y for both smooth and rough
      reward[j] <- V(g)$reward[observations[j]]
    }
    reward
  }, mc.preschedule = TRUE, mc.cores=cores)
  
  #put into dataFrame
  gpDF <- data.frame(trial=seq(1:26)-1,
                      meanReward = rowMeans(rewardTotal), #reward averaged over replications
                      meanSE = apply(rewardTotal, 1, FUN = function(x) sd(x)/sqrt(length(x))))
  #add model label
  gpDF$Model <- modelName
  #write to csv
  if(!is.null(outputfile)){
    write.csv(gpDF, outputfile)
  }
  return(gpDF)
}

gpPars <- read.csv('rationalModels/parameters/DF.csv')
gpDF <- gpRationalModel(reps, modelName = 'GP-UCB', outputfile = "rationalModels/GPUCB.csv", parameters = gpPars)

#############################################################################################################################
# dNN Model
# to read the previous simulation from disk:
# dNNDF <- read.csv("rationalModels/dNN.csv")
#############################################################################################################################

dNNRationalModel <- function(replications, modelName, outputfile, parameters, stickiness=T){
  #run for smooth environments
  rewardTotal <- mcmapply(1:replications, FUN=function(x){
    #sample parameters
    params <- sample_n(parameters,1)
    d <- params$d
    tau <- params$tau
    omega <- params$omega
    #tau <- median(smoothPars$tau)
    #randomly choose environment
    envNum <- sample(1:40,1) 
    observations <- rep(NA,26)
    #1st trial is random
    location <-observations[1]<- sample(1:64,1)#each location as an integer from 1:25; first value is random, so treat both environments as replications of the same actions
    #Observations
    g <- graphList[[envNum]]
    distanceMatrix <- distanceTables[[envNum]]
    #reward
    reward<- rep(NA, 26)
    reward[1] <- V(g)$reward[location] #first sample is random
    for (j in 2:26){ #after that, loop through remaining trials and make decisions based on BMT-UCB preditions
      #update posterior predictions
      post <- dnnChoice(distanceMatrix,X=observations[1:(j-1)], Y=(reward[1:(j-1)]- 50)/100, dist = d, mu_0=0)
      #compute acquisition function evaluation
      utilityVec <- post$mu #Mean only model since no uncertainty judgments
      if (stickiness==T){ #add stickiness
        utilityVec[location]<- utilityVec[location] + omega
      }
      #to prevent overflow, subtract max of q(x) vector 
      utilityVec <- utilityVec - max(utilityVec)
      #compute softmax choice probabilities
      p <- exp(utilityVec/tau)
      p <- p/sum(p)
      #Sample next choice
      observations[j] <-  sample(1:64,1, prob=p, replace=TRUE)
      #update reward, X1, X2, and Y for both smooth and rough
      reward[j] <- V(g)$reward[observations[j]]
    }
    reward
  }, mc.preschedule = TRUE, mc.cores=cores)
  
  #put into dataFrame
  dNNDF <- data.frame(trial=seq(1:26)-1,
                     meanReward = rowMeans(rewardTotal), #reward averaged over replications
                     meanSE = apply(rewardTotal, 1, FUN = function(x) sd(x)/sqrt(length(x))))
  #add model label
  dNNDF$Model <- modelName
  #write to csv
  if(!is.null(outputfile)){
    write.csv(dNNDF, outputfile)
  }
  return(dNNDF)
}

dNNPars <- read.csv('rationalModels/parameters/dNN.csv')
dNNDF <- dNNRationalModel(reps, modelName = 'dNN', outputfile = "rationalModels/dNN.csv", parameters = dNNPars)



#############################################################################################################################
# kNN Model
# to read the previous simulation from disk:
# kNNDF <- read.csv("rationalModels/kNN.csv")
#############################################################################################################################

kNNRationalModel <- function(replications, modelName, outputfile, parameters, stickiness=T){
  #run for smooth environments
  rewardTotal <- mcmapply(1:replications, FUN=function(x){
    #sample parameters
    params <- sample_n(parameters,1)
    k <- params$k
    tau <- params$tau
    omega <- params$omega
    #tau <- median(smoothPars$tau)
    #randomly choose environment
    envNum <- sample(1:40,1) 
    observations <- rep(NA,26)
    #1st trial is random
    location <-observations[1]<- sample(1:64,1)#each location as an integer from 1:25; first value is random, so treat both environments as replications of the same actions
    #Observations
    g <- graphList[[envNum]]
    distanceMatrix <- distanceTables[[envNum]]
    #reward
    reward<- rep(NA, 26)
    reward[1] <- V(g)$reward[location] #first sample is random
    for (j in 2:26){ #after that, loop through remaining trials and make decisions based on BMT-UCB preditions
      #update posterior predictions
      post <- knnChoice(distanceMatrix,X=observations[1:(j-1)], Y=(reward[1:(j-1)]- 50)/100, kValue = k, mu_0=0)
      #compute acquisition function evaluation
      utilityVec <- post$mu #Mean only model since no uncertainty judgments
      if (stickiness==T){ #add stickiness
        utilityVec[location]<- utilityVec[location] + omega
      }
      #to prevent overflow, subtract max of q(x) vector 
      utilityVec <- utilityVec - max(utilityVec)
      #compute softmax choice probabilities
      p <- exp(utilityVec/tau)
      p <- p/sum(p)
      #Sample next choice
      observations[j] <-  sample(1:64,1, prob=p, replace=TRUE)
      #update reward, X1, X2, and Y for both smooth and rough
      reward[j] <- V(g)$reward[observations[j]]
    }
    reward
  }, mc.preschedule = TRUE, mc.cores=cores)
  
  #put into dataFrame
  kNNDF <- data.frame(trial=seq(1:26)-1,
                      meanReward = rowMeans(rewardTotal), #reward averaged over replications
                      meanSE = apply(rewardTotal, 1, FUN = function(x) sd(x)/sqrt(length(x))))
  #add model label
  kNNDF$Model <- modelName
  #write to csv
  if(!is.null(outputfile)){
    write.csv(kNNDF, outputfile)
  }
  return(kNNDF)
}

kNNPars <- read.csv('rationalModels/parameters/kNN.csv')
kNNDF <- kNNRationalModel(reps, modelName = 'kNN', outputfile = "rationalModels/kNN.csv", parameters = kNNPars)

#############################################################################################################################
# # PLOTS
# # Load all models
randomDF <- read.csv("rationalModels/random.csv")
bmtDF <- read.csv("rationalModels/BMTUCB.csv")
gpDF <- read.csv("rationalModels/GPUCB.csv")
knnDF <- read.csv("rationalModels/kNN.csv")
dnnDF <- read.csv("rationalModels/dNN.csv")

#############################################################################################################################

#Join all DFs together 
rationalDF <- rbind(randomDF, bmtDF, gpDF, knnDF, dnnDF)
rationalDF <- rationalDF[ , !(names(rationalDF) %in% c("X"))]

#add human data
d <- dataImportBandit('../data/banditTask/networkBandit.csv', delim=',')
d <- subset(d, round<10) #remove bonus round

dplot<-ddply(d,~trial,summarise,meanReward=mean(rewardObjective), meanSE=sd(rewardObjective)/sqrt(length(rewardObjective)))
dplot$Model <- rep("Human", nrow(dplot))
#Join human with rational models
rationalDF <- rbind(rationalDF, dplot)

summary(rationalDF)
#rationalDF$trial <- rationalDF$trial - 1

#reorder factor levels
rationalDF$Model <- factor(rationalDF$Model, levels = c("Human", "GP-UCB", "BMT-UCB" ,  'dNN', 'kNN',"Random" ))
levels(rationalDF$Model) <- c("Human", "GP", "BMT" , 'dNN', 'kNN',"Rand")


#Plot of mean Reward
p1<- ggplot(rationalDF, aes(x=trial, y=meanReward, col=Model, shape = Model))+
  stat_summary(fun.y = mean, geom='line', size = 0.8)+
  #stat_summary(fun.y=mean, geom='point', size = 2)+
  geom_point(size=2) +
  #geom_line(size=.8) +
  #geom_errorbar(aes(ymin=meanReward-meanSE, ymax=meanReward+meanSE), width=0.1, size=1) +
  ylab("Average Reward")+xlab("Trial")+
  theme_classic()+
  scale_shape_manual(values = c(7, 16,17,15,4,3,32),  name='') +
  #scale_color_manual(values=c("black",  "#F0E442", "#E69F00", "#009E73", "#56B4E9", "#fb9a99"))+
  #scale_color_manual(values=c("black","#7570B3", "#1B9E77", "#fb9a99"))+
  scale_color_manual(values=c("#CC79A7",  "#009E73", "#F0E442","#e41a1c", "#56B4E9","black"), name='')+
  #scale_linetype_manual(values = c('solid', 'dotted', 'solid', 'dotted', 'solid', 'solid'))+
  #scale_color_brewer(palette="Paired", direction=1)+
  #coord_cartesian(ylim=c(45,85))+
  guides(color=guide_legend(ncol=2), shape=guide_legend(ncol=2)) +
  scale_x_continuous(breaks = round(seq(0,20, by = 5),1))+
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.05), legend.key=element_blank(),legend.background=element_blank(),
        strip.background=element_blank(), text=element_text(size=12,  family="sans"), legend.text=element_text(size=10), legend.key.size =  unit(0.15, "in"))
p1

