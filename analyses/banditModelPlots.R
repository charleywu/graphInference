#Interpret modeling results
#Charley Wu, March 2019

#house keeping
rm(list=ls())

#load packages
packages <- c("coin", 'ggbeeswarm', 'pander', 'lsr','cowplot', 'BayesFactor', 'readr', 'plyr', 'jsonlite', 'ggplot2', "ggExtra", 'gridExtra', 'reshape2', 'stargazer', 'coefplot', "grid", 'corrplot', 'ggsignif')
lapply(packages, require, character.only = TRUE)

source('utilities.R')
source('exportImportGraph.R')
source('statisticalTests.R')

theme_set(theme_bw(base_size=16))# use the b&w theme
#############################################################################################################################
# IMPORT PARTICIPANT DATA 
#############################################################################################################################

#Participant data
data <- dataImportBandit('../data/banditTask/networkBandit.csv', delim=',')
data <- subset(data, round<10) #remove bonus round
uids <- unique(data$id)


#Wrapper for brm models such that it saves the full model the first time it is run, otherwise it loads it from disk
run_model <- function(expr, modelName, path='brmsModels', reuse = TRUE) {
  path <- paste0(path,'/', modelName, ".brm")
  if (reuse) {
    fit <- suppressWarnings(try(readRDS(path), silent = TRUE))
  }
  if (is(fit, "try-error")) {
    fit <- eval(expr)
    saveRDS(fit, file = path)
  }
  fit
}
#############################################################################################################################
# IMPORT MODELING RESULTS
#############################################################################################################################

importModelResults <- function(dataFolder, kernels, acqFuncs){
  #initialize data frames
  modelFit <- data.frame() 
  paramEstimates <- data.frame()
  choicePredictions <- data.frame(participant=numeric(), round=numeric(), trial = numeric(), choice=numeric(),  kernel=numeric(), acq=numeric(), nLL=numeric())
  #loop through data
  for (k in kernels){
    for (a in acqFuncs){
      for (i in uids){ #subjects
        filename <- paste0(dataFolder,"/", k, a, i, ".Rdata") #read output file
        if (file.exists(filename)){
          #debug: 
          #filename <- 'modelResults/networkBandit/DFUCB4.Rdata'
          #k <- 'kNN'
          #a <- 'GM'
          #i <- 54
          datalist<-get(load(filename))
          dp<- datalist[[1]] #MLE
          #Add missing columns
          colnames <- c("loo", "nLL", "tau", "beta", "lambda", "errorVariance", "alpha", "d", "k",  "omega")
          add <-colnames[!colnames%in%names(dp)]
          if(length(add)!=0) dp[add] <- NA
          #names(dp) <- c("loo", "nLL", "tau", "beta", "lambda", "errorVariance", "alpha") #DEBUG
          smlist <- datalist[[2]] #Softmax prediction surface
          choiceHistory <- datalist[[3]] #choices
          #demographics
          dummy<-subset(data, id==i) #subset participant in the dataframe
          environment <- dummy$environment[1]
          participant<-dummy$id[1]  #replicate ID
          kernel<-k
          acq<-a
          #Total out of sample nLL
          nLL <- sum(dp$nLL)
          randomModelLL <- -log(1/64)*225
          R2 <- 1 - (nLL/randomModelLL)
          #save modelFit
          dadd<-data.frame(participant=participant, nLL=nLL, kernel=kernel, acq = acq, R2 =R2, kError = median(dp$errorVariance), lambda=median(dp$lambda), beta = median(dp$beta), tau=median(dp$tau), alpha = median(dp$alpha), d = median(dp$d), k=median(dp$k), omega = median(dp$omega) )
          modelFit<-rbind(modelFit, dadd)
          #Save all LOO parameter estimates for each subject
          dp$participant <- participant
          dp$R2 <- 1 - (dp$nLL / (-log(1/64)*200))
          dp$acq <- a
          dp$kernel <- k
          paramEstimates<-rbind(paramEstimates, dp)
          #Save log loss for individual choices
          individualChoices <- sapply(1:9, FUN= function(r){  #loop through rounds
            roundM <- smlist[[r]] #choice matrix for round r
            choices <- cbind(1:25, choiceHistory[r,1:25]) #join into two column matrix: trial number, choice selection
            nLL <- - log(roundM[choices]) #compute log loss for each individual choice
          } ) #returns trial x round matrix of log loss values for each individual choice
          choiceVector <- as.vector(individualChoices) #unravel into vector
          roundVec <- rep(seq(1,9), each=25)
          trialVec <- rep(seq(1,25), 9)
          dummyDF <- data.frame(participant=participant, round= roundVec, trial = trialVec, choice = as.vector(choiceHistory[,1:25]),  kernel=kernel, acq = acq, nLL = choiceVector)
          choicePredictions <- rbind(choicePredictions, dummyDF)
          }}}}
  return(list(modelFit, paramEstimates, choicePredictions))
}

#############################################################################################################################
# COMPILE MODELING RESULTS FROM POTENTIALLY MULTIPLE SOURCES
# Load all models so far:
#modelFit <- read.csv('modelResults/modelFit.csv')
#paramEstimates <- read.csv('modelResults/paramEstimates.csv')
#############################################################################################################################
modelList <- c("BMT",  'DF', 'kNN', 'dNN', 'Sticky')
acqusitionFuncs<-  c("UCB", 'GM')

modelResults <- importModelResults('modelResults/networkBandit/', modelList, acqusitionFuncs)

#separate into overall per participant and individual cross validation blocks
modelFit <- modelResults[[1]]
paramEstimates <- modelResults[[2]]
choicePredictions <- modelResults[[3]]

#reorder acqisition function levels and add "ModelName" to identify unique models
modelFit$acq <-factor(modelFit$acq)
modelFit$kernel <- factor(modelFit$kernel)
#levels(modelFit$kernel) <- c("BMT", "RBF")
modelFit$ModelName <- paste(modelFit$kernel,modelFit$acq, sep = "-")
modelFit$ModelName <- factor(modelFit$ModelName)
unique(modelFit$ModelName)

paramEstimates$kernel <- factor(paramEstimates$kernel)
#levels(paramEstimates$kernel) <- c("BMT", "RBF")

#add participant data to modelFit
meanDF <- ddply(data, .(id), plyr::summarize, meanScore = mean(rewardObjective))
modelFit <- merge(modelFit, meanDF, by.x=c("participant"), by.y=c("id"))

#add participant data to paramEstimates
paramEstimates <- merge(paramEstimates, meanDF, by.x=c("participant"), by.y=c("id"))

#SAVE MODEL RESULTS TO DISK
#write.csv(modelFit,'modelResults/modelFit.csv')
#write.csv(paramEstimates,'modelResults/paramEstimates.csv')

#Save nLL for PXP
write.table(cbind(subset(modelFit, ModelName == 'DF-UCB')$nLL, subset(modelFit,  ModelName == 'BMT-UCB')$nLL, subset(modelFit,  ModelName == 'dNN-GM')$nLL, subset(modelFit,  ModelName == 'kNN-GM')$nLL,  subset(modelFit,  ModelName == 'Sticky-GM')$nLL), 'modelResults/Exp2diffevidence.csv', row.names = FALSE, sep=',', col.names=FALSE)

#############################################################################################################################
# SAVE PARAMETER ESTIMATES FOR RATIONAL MODELS USED TO COMPUTE LEARNING CURVES
#############################################################################################################################
#BMTpars <- subset(paramEstimates, kernel=="BMT" & acq=="UCB")
#write.csv(BMTpars[,c("participant", "loo", "errorVariance", "beta", "tau", "omega")], file = "rationalModels/parameters/BMT.csv", row.names = FALSE)

#DFpars <- subset(paramEstimates, kernel=="DF" & acq=="UCB")
#write.csv(DFpars[,c("participant", "loo", "alpha", "beta", "tau", "omega")], file = "rationalModels/parameters/DF.csv", row.names = FALSE)

#kNNpars <- subset(paramEstimates, kernel=="kNN" & acq=="GM")
#write.csv(kNNpars[,c("participant", "loo", "k", "tau", "omega")], file = "rationalModels/parameters/kNN.csv", row.names = FALSE)

#dNNpars <- subset(paramEstimates, kernel=="dNN" & acq=="GM")
#write.csv(dNNpars[,c("participant", "loo", "d", "tau", "omega")], file = "rationalModels/parameters/dNN.csv", row.names = FALSE)
#############################################################################################################################
# Summary of Modeling Results
#############################################################################################################################
#ALL MODELS
#Number of participants best described
models <-  rep(0, length(levels(modelFit$ModelName)))
names(models) <- levels(modelFit$ModelName)
uids <- unique(subset(modelFit)$participant)
for (pid in uids){
  subDF <- subset(modelFit, participant==pid)
  best <- subDF$ModelName[which(subDF$R2==max(subDF$R2))]
  models[best] <- models[best] + 1
}
#add best described to modelFit 
bestDescribed <- paste0('bestDescribed')
modelFit[bestDescribed] <- NA
for (mod in names(models)){
  modelFit[modelFit$ModelName == mod,bestDescribed] <- models[mod]
}
#Summary of results

for (k_i in levels(modelFit$kernel)){
  for (a in levels(modelFit$acq)){
    if (paste(k_i,a,sep="-") %in% names(models)){
      sub <- subset(modelFit, acq == a & kernel == k_i) #sub of model fit
      psub <- subset(paramEstimates, acq == a & kernel==k_i)
      randomModel <- -log(1/64) * 225
      r_2 <- 1 - (mean(sub$nLL)/randomModel)
      bestDescribed <- models[paste(k_i,a,sep="-")]
      cat( k_i,"-",a, " (n=",nrow(sub),")", "; R2 = ", round(r_2, digits=2),"; Best Described = ", bestDescribed , "; Alpha = ", round(median(psub$alpha), digits=2),
           "; Lambda = ", round(median(psub$lambda), digits=2),"; Beta = ", round(median(psub$beta), digits=2), "; kError = ", round(median(psub$errorVariance), digits=2), 
           "; k = ", round(median(psub$k), digits=2), "; d = ", round(median(psub$d), digits=2), "; gamma = ", round(median(psub$gamma), digits=2), 
           "; Tau = ", round(median(psub$tau), digits=2),  "; omega = ", round(median(psub$omega), digits=2),"\n",  sep="")
    }
  }
}
cat("\n")

  
#############################################################################################################################
# Summary of Modeling Results
#############################################################################################################################
#modelFit <- subset(modelFit, kernel != 'Sticky')
modelFit$kernel <- factor(modelFit$kernel,levels = c("DF","BMT","dNN", "kNN", "Sticky"))
levels(modelFit$kernel) <-c("GP","BMT", "dNN", "kNN", "Sticky")
plotDF <- modelFit
plotDF$kernel <- factor(plotDF$kernel)


se<-function(x){sd(x)/sqrt(length(x))}
mtmDF <- ddply(modelFit, ~participant+kernel, summarise, newR2 =mean(R2), se=se(R2))
boxplotDF <- ddply(plotDF, ~kernel, summarise,d_ymin = max(min(R2), quantile(R2, 0.25) - 1.5 * IQR(R2)), d_ymax = min(max(R2), quantile(R2, 0.75) + 1.5 * IQR(R2)),
                   d_lower = quantile(R2, 0.25),  d_middle = median(R2), d_upper = quantile(R2, 0.75),
                   mu=mean(R2))

p1 <- ggplot(boxplotDF) +
  #stat_summary(fun.y = mean, geom = "bar", position = "dodge", color='black') + 
  geom_boxplot(aes(x =0-0.2, ymin = d_lower, ymax = d_upper, lower = d_lower, middle = d_middle, upper = d_upper, width = 2 * 0.2, fill = kernel), stat = "identity", color='black') +
  #whiskers
  geom_segment(aes(x = 0, y = d_ymin, xend = 0, yend = d_ymax)) +
  geom_segment(aes(x = 0 - 0.1,  y = d_ymax - 0.0025, xend = 0,  yend = d_ymax -0.0025)) + #top
  geom_segment(aes(x = 0 - 0.1, y = d_ymin +0.0025, xend = 0, yend = d_ymin + 0.0025)) + #bottom
  geom_point(aes(x = 0-0.2, y = mu), size = 2, shape=23, fill='white') +
  geom_jitter(data=mtmDF, aes(x = 0+.2,  y = newR2,  color = kernel), 
              width = 0.2 - 0.25 * 0.2, height = 0, size=1, alpha = 0.7)+
  #geom_jitter(data = mainTextModels, aes(x=acq, y = R2, fill= kernel),  color='grey', size = 1, shape=21, alpha=0.2, position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.2))+
  #geom_errorbar(aes(ymin=newR2 - se, ymax=newR2 + se),color='black', width = .4, position=position_dodge((width=0.9))) +
  #scale_fill_manual(values = c( "#009E73", "#E69F00", "#56B4E9"), name ="")+
  #scale_color_manual(values = c( "#009E73", "#E69F00", "#56B4E9"), name="") +
  scale_fill_manual(values = c( "#009E73", "#F0E442","#E69F00", "#56B4E9","#e31a1c"), name ="")+
  scale_color_manual(values = c( "#009E73", "#F0E442", "#E69F00", "#56B4E9","#e31a1c"), name="") +
  #scale_color_brewer(palette = 'Set1', direction = -1, name = '')+
  #scale_fill_brewer(palette = 'Set1', direction = -1, name = '')+
  ylab("Predictive Accuracy")+ 
  xlab('Task')+
  facet_grid(~kernel, scales='free_x')+
  #coord_cartesian(ylim=c(-0.22, 0.8))+
  theme_classic()+
  theme(text = element_text(size=12,  family="sans"),
        strip.background=element_blank(),
        #strip.text = element_blank(),
        legend.key=element_rect(color=NA),
        panel.spacing.x=unit(0.2, "lines"),
        panel.spacing.y=unit(1, "lines"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position='none')+
  guides(color=FALSE, shape=FALSE)
p1


#t-tests
ttestPretty(subset(modelFit, kernel=='GP')$R2-subset(modelFit, kernel=='dNN')$R2, mu=0)
ttestPretty(subset(modelFit, kernel=='GP')$R2-subset(modelFit, kernel=='BMT')$R2, mu=0, maxBF=10000)
ttestPretty(subset(modelFit, kernel=='GP')$R2-subset(modelFit, kernel=='kNN')$R2, mu=0, maxBF=10000)

ttestPretty(subset(modelFit, kernel=='dNN')$R2-subset(modelFit, kernel=='BMT')$R2, mu=0, maxBF=10000)
ttestPretty(subset(modelFit, kernel=='kNN')$R2-subset(modelFit, kernel=='BMT')$R2, mu=0, maxBF=10000)
ttestPretty(subset(modelFit, kernel=='dNN')$R2-subset(modelFit, kernel=='kNN')$R2, mu=0, maxBF=10000)


#PXP plot (calculated from modelfit.m)
PXP <- as.numeric(read.csv('modelResults/Exp2PXP.csv', header=F))
pxpDF <- data.frame(model= c('GP', 'BMT', 'dNN', 'kNN', 'Sticky'),pxp= PXP) #build dataframe with output
#Reorder levels
pxpDF$model <- factor(pxpDF$model , levels=c('GP', 'BMT', 'dNN', 'kNN', 'Sticky'))

p1b <- ggplot(pxpDF, aes(x = model, y = pxp, fill = model))+
  geom_bar(stat="identity",color='black')+
  ylab('pxp')+
  scale_fill_manual(values = c( "#009E73", "#F0E442","#E69F00", "#56B4E9","#e41a1c"), name ="")+
  xlab('')+
  scale_y_continuous(labels = percent)+
  expand_limits(y=0)+
  theme_classic()+
  ylim(c(0,1))+
  #scale_fill_manual(values = c(gpCol,"#a6761d", bmtCol), name='')+
  #scale_fill_manual(values = c(gpCol,bmtCol), name='')+
  theme(legend.position="none", strip.background=element_blank(), legend.key=element_rect(color=NA))
p1b


#Diagnosticity between GP and BMT
subset(modelFit, kernel %in% c('GP', 'BMT')) %>% group_by(participant,kernel, meanScore) %>%summarize(diff=diff(nLL))
modelDiff <- ddply(subset(modelFit, kernel %in% c('GP', 'BMT')), ~participant, summarize, meanScore = mean(meanScore))
modelDiff$diff <- subset(modelFit, kernel=='BMT')$nLL - subset(modelFit, kernel=='GP')$nLL
modelDiff$color <- ifelse(modelDiff$diff>0,"GP", "BMT" )
modelDiff$color <- factor(modelDiff$color, levels = c('GP', 'BMT'))

corTestPretty(modelDiff$meanScore, modelDiff$diff)
pModelDiffGP_BMT <- ggplot(modelDiff, aes(x = meanScore, y = diff ))+
  geom_point(aes(color = color), alpha = .9)+
  geom_hline(yintercept = 0, color ='black', linetype = 'dashed')+
  geom_smooth(method='lm', alpha = 0.2, color = 'black')+
  theme_classic() +
  xlab('Mean Bandit Score')+
  ylab(expression(nLL[BMT] - nLL[GP]))+
  scale_fill_manual(values = c( "#009E73", "#F0E442","#E69F00", "#56B4E9","#e31a1c"), name ="Better Model")+
  scale_color_manual(values = c( "#009E73", "#F0E442", "#E69F00", "#56B4E9","#e31a1c"), name ="Better Model") +
  annotate("text", x = 70, y = 135, label = "paste(italic(r) , \" = .21\", \", \", italic(BF)[10] == 1.8 )", parse = TRUE)+
  theme(legend.position=c(0,0), legend.justification = c(0,0),strip.background=element_blank(), legend.key=element_blank(), legend.background=element_blank())
pModelDiffGP_BMT

#Diagnosticity between GP and dNN
subset(modelFit, kernel %in% c('GP', 'dNN')) %>% group_by(participant,kernel, meanScore) %>%summarize(diff=diff(nLL))
modelDiff <- ddply(subset(modelFit, kernel %in% c('GP', 'dNN')), ~participant, summarize, meanScore = mean(meanScore))
modelDiff$diff <- subset(modelFit, kernel=='dNN')$nLL - subset(modelFit, kernel=='GP')$nLL
modelDiff$color <- ifelse(modelDiff$diff>0,"GP", "dNN" )
modelDiff$color <- factor(modelDiff$color, levels = c('GP', 'dNN'))

corTestPretty(modelDiff$meanScore, modelDiff$diff)
pModelDiffGP_dNN <- ggplot(modelDiff, aes(x = meanScore, y = diff ))+
  geom_point(aes(color = color), alpha = .9)+
  geom_hline(yintercept = 0, color ='black', linetype = 'dashed')+
  geom_smooth(method='lm', alpha = 0.2, color = 'black')+
  theme_classic() +
  xlab('Mean Bandit Score')+
  ylab(expression(nLL[dNN] - nLL[GP]))+
  scale_fill_manual(values = c( "#009E73", "#E69F00", "#56B4E9","#e31a1c"), name ="Better Model")+
  scale_color_manual(values = c( "#009E73", "#E69F00", "#56B4E9","#e31a1c"), name ="Better Model") +
  annotate("text", x = 70, y = 100, label = "paste(italic(r) , \" = -.07\", \", \", italic(BF)[10] == 0.29 )", parse = TRUE)+
  theme(legend.position=c(0,0), legend.justification = c(0,0),strip.background=element_blank(), legend.key=element_blank(), legend.background=element_blank())
pModelDiffGP_dNN


#Diagnosticity between GP and kNN
subset(modelFit, kernel %in% c('GP', 'kNN')) %>% group_by(participant,kernel, meanScore) %>%summarize(diff=diff(nLL))
modelDiff <- ddply(subset(modelFit, kernel %in% c('GP', 'kNN')), ~participant, summarize, meanScore = mean(meanScore))
modelDiff$diff <- subset(modelFit, kernel=='kNN')$nLL - subset(modelFit, kernel=='GP')$nLL
modelDiff$color <- ifelse(modelDiff$diff>0,"GP", "kNN" )
modelDiff$color <- factor(modelDiff$color, levels = c('GP', 'kNN'))

corTestPretty(modelDiff$meanScore, modelDiff$diff)
pModelDiffGP_kNN <- ggplot(modelDiff, aes(x = meanScore, y = diff ))+
  geom_point(aes(color = color), alpha = .9)+
  geom_hline(yintercept = 0, color ='black', linetype = 'dashed')+
  geom_smooth(method='lm', alpha = 0.2, color = 'black')+
  theme_classic() +
  xlab('Mean Bandit Score')+
  ylab(expression(nLL[kNN] - nLL[GP]))+
  scale_fill_manual(values = c( "#009E73",  "#56B4E9","#e31a1c"), name ="Better Model")+
  scale_color_manual(values = c( "#009E73",  "#56B4E9","#e31a1c"), name ="Better Model") +
  annotate("text", x = 70, y = 135, label = "paste(italic(r) , \" = -.03\", \", \", italic(BF)[10] == 0.24 )", parse = TRUE)+
  theme(legend.position=c(0,0), legend.justification = c(0,0),strip.background=element_blank(), legend.key=element_blank(), legend.background=element_blank())
pModelDiffGP_kNN


#Diagnosticity between GP and Sticky
subset(modelFit, kernel %in% c('GP', 'Sticky')) %>% group_by(participant,kernel, meanScore) %>%summarize(diff=diff(nLL))
modelDiff <- ddply(subset(modelFit, kernel %in% c('GP', 'Sticky')), ~participant, summarize, meanScore = mean(meanScore))
modelDiff$diff <- subset(modelFit, kernel=='Sticky')$nLL - subset(modelFit, kernel=='GP')$nLL
modelDiff$color <- ifelse(modelDiff$diff>0,"GP", "Sticky" )
modelDiff$color <- factor(modelDiff$color, levels = c('GP', 'Sticky'))

corTestPretty(modelDiff$meanScore, modelDiff$diff)
pModelDiffGP_Sticky <- ggplot(modelDiff, aes(x = meanScore, y = diff ))+
  geom_point(aes(color = color), alpha = .9)+
  geom_hline(yintercept = 0, color ='black', linetype = 'dashed')+
  geom_smooth(method='lm', alpha = 0.2, color = 'black')+
  theme_classic() +
  xlab('Mean Bandit Score')+
  ylab(expression(nLL[Sticky] - nLL[GP]))+
  scale_fill_manual(values = c( "#009E73","#e31a1c"), name ="Better Model")+
  scale_color_manual(values = c( "#009E73", "#e31a1c"), name ="Better Model") +
  annotate("text", x = 60, y = 200, label = "paste(italic(r) , \" = .41\", \", \", italic(BF)[10] > 100 )", parse = TRUE)+
  theme(legend.position=c(0,0), legend.justification = c(0,0),strip.background=element_blank(), legend.key=element_blank(), legend.background=element_blank())
pModelDiffGP_Sticky


###########################################################################################################################
#Learning Curves: Computed in modelSimulations.R
###########################################################################################################################
randomDF <- read.csv("rationalModels/random.csv")
bmtDF <- read.csv("rationalModels/BMTUCB.csv")
gpDF <- read.csv("rationalModels/GPUCB.csv")
knnDF <- read.csv("rationalModels/kNN.csv")
dnnDF <- read.csv("rationalModels/dNN.csv")


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

#reorder factor levels

rationalDF$Model <- factor(rationalDF$Model, levels = c("Random",  "GP-UCB", "BMT-UCB" , 'dNN', 'kNN', "Human"))
levels(rationalDF$Model) <- c("Random",  "GP", "BMT" , 'dNN', 'kNN',"Human")


#Plot of mean Reward
pLearningCurve<- ggplot(rationalDF, aes(x=trial, y=meanReward, col=Model, shape = Model))+
  stat_summary(fun.y = mean, geom='line', size = 0.8)+
  #stat_summary(fun.y=mean, geom='point', size = 2)+
  geom_point(size=2) +
  #geom_line(size=.8) +
  #geom_errorbar(aes(ymin=meanReward-meanSE, ymax=meanReward+meanSE), width=0.1, size=1) +
  ylab("Average Reward")+xlab("Trial")+
  theme_classic()+
  scale_shape_manual(values = c(32, 16,17,15,4,3,7), name='') +
  #scale_color_manual(values=c("black",  "#F0E442", "#E69F00", "#009E73", "#56B4E9", "#fb9a99"))+
  #scale_color_manual(values=c("black","#7570B3", "#1B9E77", "#fb9a99"))+
  scale_color_manual(values=c( "black", "#009E73", "#F0E442","#E69F00", "#56B4E9", "#CC79A7"), name='')+
  #scale_linetype_manual(values = c('solid', 'dotted', 'solid', 'dotted', 'solid', 'solid'))+
  #scale_color_brewer(palette="Paired", direction=1)+
  #coord_cartesian(ylim=c(45,85))+
  guides(color=guide_legend(ncol=2), shape=guide_legend(ncol=2)) +
  scale_x_continuous(breaks = round(seq(0,20, by = 5),1))+
  theme(legend.justification = c(1, 0), legend.position = c(1, 0.05), legend.key=element_blank(),legend.background=element_blank(),
        strip.background=element_blank(), text=element_text(size=12,  family="sans"), legend.text=element_text(size=10), legend.key.size =  unit(0.15, "in"))
pLearningCurve


###########################################################################################################################
#Parameter Estimates
###########################################################################################################################

#Set the upper limit based on Tukey's outlier removal
upperlimit <- function(k, a, params){ #k is kernel, a is acquisition function, params is a vector of parameter names
  upperquartiles <- sapply(params, FUN=function(i) quantile(modelFit[modelFit$kernel %in% k & modelFit$acq==a, i], probs=c(.25, .75))[2]) #upper quartiles
  H <- sapply(params, FUN=function(i) 1.5 * IQR(modelFit[modelFit$kernel %in% k & modelFit$acq==a, i])) #upper quartiles
  upperLimit <- max(upperquartiles + H) #use the larger 1.5 * IQR from the upper quartile to set axis limits
  return(upperLimit)
}

colorpicker <- function(paramList){
  allParams <- c("alpha", "gamma", "kError", "d", "k",   "beta", "tau", "omega")
  dark2 <- c('#e7298a','#66a61e','#e6ab02','#a6761d','#666666','#1b9e77','#d95f02','#7570b3')
  colors <- dark2[which(allParams %in% paramList )]
  return(colors)
}

#DIFFUSION KERNEL
#Melt to separate out parameters
mDF <- melt(modelFit, id.vars = c("acq", "kernel", "R2", "participant"), measure.vars = c("alpha", "beta", "tau", 'omega' ))
#Select which acq and which kernel to display

k<- "GP"
a<- "UCB"

p2a<- ggplot(subset(mDF, kernel %in% k & acq==a), aes(y=value, x=0, color= variable)) +
  geom_quasirandom(aes(fill=variable), alpha=.5, width = 0.5) +
  geom_boxplot(aes(color=variable),color='black', fill=NA,width=0.2, outlier.shape=NA, lwd=0.5, position = 'dodge') +
  stat_summary(color='black', fill='white',fun.y=mean, geom="point", shape=23, size=2, stroke = 1) + #mean
  xlab("") +
  ylab("Estimate")+ 
  scale_fill_manual(values = colorpicker(c("alpha", "beta", "tau", 'omega' )), name = '' ) +
  scale_color_manual(values = colorpicker(c("alpha", "beta", "tau", 'omega' )), name = '' ) +
  scale_y_continuous(trans="log10") +
  annotation_logticks(base=10, sides='l')+
  #coord_cartesian(ylim=c(0,upperlimit(k,a,c("alpha", "beta", "tau", "omega" ))))+
  facet_grid(variable~., scales='free', labeller = label_parsed)+
  #facet_grid(~variable, labeller = label_parsed)+
  #coord_cartesian(ylim=c(0,2))+
  ggtitle("GP") +
  theme_classic()+
  theme(text = element_text(size=12,  family="sans"),strip.background=element_blank(), legend.key=element_rect(color=NA), legend.position="none", panel.spacing.x=unit(0.2, "lines"),
        panel.spacing.y=unit(1, "lines"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text.y = element_text(angle = 0, size = 14),
        axis.ticks.x=element_blank()) 
#scale_x_discrete("",labels=c("")) + 
p2a
ttestPretty(subset(modelFit, kernel =='GP')$alpha, mu=2, maxBF=Inf)

ttestPretty(subset(mDF, kernel %in% k & acq==a & variable=='beta')$value, mu=exp(-5), maxBF=Inf)
corTestPretty(subset(mDF, kernel %in% k & acq==a & variable=='beta')$value, subset(mDF, kernel %in% k & acq==a & variable=='tau')$value, method='kendall')
corTestPretty(subset(mDF, kernel %in% k & acq==a & variable=='alpha')$value, subset(mDF, kernel %in% k & acq==a & variable=='beta')$value, method='kendall')


#BMT
#Melt to separate out parameters
correctedDF <- modelFit
correctedDF$kError <- correctedDF$kError^2 #kError is actually a stddev, so here we make it a variance 
mDF <- melt(correctedDF, id.vars = c("acq", "kernel", "R2", "participant"), measure.vars = c("kError", "beta", "tau", "omega"))

#Select which acq and which kernel to display
k<- "BMT"
a<- "UCB"
options(digits=4)
scaleFUN <- function(x) sprintf("%.2f", x)
levels(mDF$variable) <- c(expression("sigma[epsilon]^2"), "beta", "tau", "omega")
p2c<-  ggplot(subset(mDF, kernel %in% k & acq==a), aes(y=value, x=0, color= variable)) +
  geom_quasirandom(aes(fill=variable), alpha=.5, width = 0.5) +
  geom_boxplot(aes(color=variable),color='black', fill=NA,width=0.2, outlier.shape=NA, lwd=0.5, position = 'dodge') +
  stat_summary(color='black', fill='white',fun.y=mean, geom="point", shape=23, size=2, stroke = 1) + #mean
  xlab("") +
  ylab("")+ 
  scale_fill_manual(values =c('#e6ab02', '#1b9e77','#d95f02','#7570b3'), name = '' ) +
  scale_color_manual(values = c('#e6ab02', '#1b9e77','#d95f02','#7570b3'), name = '' ) +
  scale_y_continuous(trans="log10", labels=scaleFUN) +
  annotation_logticks(base=10, sides='l')+
  #coord_cartesian(ylim=c(0,upperlimit(k,a,c("kError", "beta", "tau", "omega"))))+
  facet_grid(variable~., scales='free', labeller = label_parsed)+
  #coord_cartesian(ylim=c(0,2))+
  ggtitle("BMT") +
  theme_classic()+
  theme(text = element_text(size=14,  family="sans"),strip.background=element_blank(), legend.key=element_rect(color=NA), legend.position="none", panel.spacing.x=unit(0.2, "lines"),
        panel.spacing.y=unit(1, "lines"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text.y = element_text(angle = 0, size = 14),
        axis.ticks.x=element_blank()) 
p2c
ttestPretty(subset(modelFit, kernel =='BMT')$beta, mu=exp(-5))

#dNN
#Melt to separate out parameters
mDF <- melt(modelFit, id.vars = c("acq", "kernel", "R2", "participant"), measure.vars = c( "d", "tau", "omega"))
#Select which acq and which kernel to display
k<- "dNN"
a<- "GM"


p2d<-  ggplot(subset(mDF, kernel %in% k & acq==a), aes(y=value, x=0, color= variable)) +
  geom_quasirandom(aes(fill=variable), alpha=.5, width = 0.5) +
  geom_boxplot(aes(color=variable),color='black', fill=NA,width=0.2, outlier.shape=NA, lwd=0.5, position = 'dodge') +
  stat_summary(color='black', fill='white',fun.y=mean, geom="point", shape=23, size=2, stroke = 1) + #mean
  xlab("") +
  ylab("")+ 
  scale_fill_manual(values = colorpicker(c("d",  "tau", 'omega' )), name = '' ) +
  scale_color_manual(values = colorpicker(c("d",  "tau", 'omega' )), name = '' ) +
  scale_y_continuous(trans="log10") +
  annotation_logticks(base=10, sides='l')+
  #coord_cartesian(ylim=c(0,upperlimit(k,a,c("kError", "beta", "tau", "omega"))))+
  facet_grid(variable~., scales='free', labeller = label_parsed)+
  #facet_grid_sc(rows = vars(variable), scales = list(y = scales_y))+
  #coord_cartesian(ylim=c(0,2))+
  ggtitle("dNN") +
  theme_classic()+
  theme(text = element_text(size=14,  family="sans"),strip.background=element_blank(), legend.key=element_rect(color=NA), legend.position="none", panel.spacing.x=unit(0.2, "lines"),
        panel.spacing.y=unit(1, "lines"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text.y = element_text(angle = 0, size = 14),
        axis.ticks.x=element_blank()) 
p2d
ttestPretty(subset(modelFit, kernel=='dNN')$tau, subset(modelFit, kernel=='DF')$tau, paired=T)
ttestPretty(subset(modelFit, kernel=='dNN')$omega, subset(modelFit, kernel=='DF')$omega, paired=T)
#kNN

#Melt to separate out parameters
mDF <- melt(modelFit, id.vars = c("acq", "kernel", "R2", "participant"), measure.vars = c( "k", "tau", "omega"))
#Select which acq and which kernel to display
k<- "kNN"
a<- "GM"

p2e<-  ggplot(subset(mDF, kernel %in% k & acq==a), aes(y=value, x=0, color= variable)) +
  geom_quasirandom(aes(fill=variable), alpha=.5, width = 0.5) +
  geom_boxplot(aes(color=variable),color='black', fill=NA,width=0.2, outlier.shape=NA, lwd=0.5, position = 'dodge') +
  stat_summary(color='black', fill='white',fun.y=mean, geom="point", shape=23, size=2, stroke = 1) + #mean
  xlab("") +
  ylab("")+ 
  scale_fill_manual(values = colorpicker(c("k",  "tau", 'omega' )), name = '' ) +
  scale_color_manual(values = colorpicker(c("k",  "tau", 'omega' )), name = '' ) +
  scale_y_continuous(trans="log10") +
  annotation_logticks(base=10, sides='l')+
  #coord_cartesian(ylim=c(0,upperlimit(k,a,c("kError", "beta", "tau", "omega"))))+
  facet_grid(variable~., scales='free', labeller = label_parsed)+
  #coord_cartesian(ylim=c(0,2))+
  ggtitle("kNN") +
  theme_classic()+
  theme(text = element_text(size=14,  family="sans"),strip.background=element_blank(), legend.key=element_rect(color=NA), legend.position="none", panel.spacing.x=unit(0.2, "lines"),
        panel.spacing.y=unit(1, "lines"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        strip.text.y = element_text(angle = 0, size = 14),
        axis.ticks.x=element_blank()) 
p2e


pParams <- cowplot::plot_grid(p2a, p2c,  p2d, p2e, nrow = 1)
pParams


mean(subset(modelFit,kernel=='kNN')$k)
median(subset(modelFit,kernel=='kNN')$k)

#Compare BMT to GP tau
ttestPretty(subset(modelFit,kernel=='BMT')$tau, subset(modelFit,kernel=='GP')$tau, paired=T)
#Compare GP stickiness to BMT stickiness
ttestPretty(subset(modelFit,kernel=='BMT')$omega, subset(modelFit,kernel=='GP')$omega, paired=T)
#random exploration of NN models
ttestPretty(subset(modelFit,kernel=='kNN')$tau, subset(modelFit,kernel=='dNN')$tau, paired=T)
#stickiness NN models
ttestPretty(subset(modelFit,kernel=='kNN')$omega, subset(modelFit,kernel=='dNN')$omega, paired=T)
###########################################################################################################################
#Performance correlation
###########################################################################################################################
#
for (model in unique(modelFit$kernel)){
  cat(model)
  modelSub<- subset(modelFit, kernel==model)
  print(corTestPretty(modelSub$R2, modelSub$meanScore))
}


#Accuracy and performance
corTestPretty(subset(modelFit, kernel=='GP')$nLL,subset(modelFit, kernel=='GP')$meanScore, maxBF=Inf )
pAccuracyR2 <- ggplot(modelFit, aes(x = meanScore, y = R2, color = kernel, fill=kernel, shape=kernel ))+
  geom_point(alpha= 0.9)+
  #geom_smooth(method='lm')+
  theme_classic() +
  xlab('Mean Bandit Score')+
  scale_shape_manual(values=c(16, 17, 15, 14, 18), name='')+
  scale_fill_manual(values = c( "#009E73", "#F0E442","#E69F00", "#56B4E9","#e31a1c"), name ="")+
  scale_color_manual(values = c( "#009E73", "#F0E442", "#E69F00", "#56B4E9","#e31a1c"), name="") +
  ylab('Predictive Accuracy')+
  #scale_color_manual(values = c(gpCol,bmtCol), name='')+
  theme(legend.position=c(0,1),legend.justification=c(0,1), strip.background=element_blank(), legend.background = element_blank())
pAccuracyR2

corTestPretty(subset(modelFit, kernel=='GP')$alpha ,subset(modelFit, kernel=='GP')$R2, maxBF=Inf )
p3a<- ggplot(subset(modelFit, kernel =='GP'), aes(x=alpha, y=R2, color=kernel, fill=kernel, shape = kernel))+
  #geom_jitter(alpha=0.4)+
  stat_smooth(method='lm', alpha = 0.1, fullrange=TRUE) +
  geom_point(alpha=0.7)+
  #geom_rug(alpha=0.5)+
  theme_classic()+
  xlim(c(0, upperlimit('GP', 'UCB', c('alpha'))))+ #Max axis limit based on the Tukey fences
  #scale_x_continuous(trans=log10_trans(), breaks = c(0.001, 0.01, 0.1, 1, 10, 100), labels=  c(0.001, 0.01, 0.1, 1, 10, 100))+
  #annotation_logticks(sides="b")+
  #facet_grid(~locality, scales='free')+
  scale_shape_manual(values=c(16, 17, 15, 16, 17), name='')+
  scale_fill_manual(values = c("#009E73", "#F0E442","#E69F00", "#56B4E9","#e31a1c"), name ="")+
  scale_color_manual(values = c("#009E73", "#F0E442","#E69F00", "#56B4E9","#e31a1c"), name="") +
  #ggtitle(bquote(alpha ~~ R^2))+
  ylab(expression(R^2))+
  xlab(expression(alpha))+
  theme(legend.position=c(1,0), legend.justification=c(1,0), strip.background=element_blank(), legend.background=element_blank())
p3a

corTestPretty(subset(modelFit, kernel=='BMT')$kError ,subset(modelFit, kernel=='BMT')$R2, maxBF=Inf )
p3b<- ggplot(subset(modelFit, kernel =='BMT'), aes(x=kError, y=R2, color=kernel, fill=kernel, shape = kernel))+
  #geom_jitter(alpha=0.4)+
  stat_smooth(method='lm', alpha = 0.1, fullrange=TRUE) +
  geom_point(alpha=0.7)+
  #geom_rug(alpha=0.5)+
  theme_classic()+
  xlim(c(0, upperlimit('BMT', 'UCB', c('kError'))))+ #Max axis limit based on the Tukey fences
  #scale_x_continuous(trans=log10_trans(), breaks = c(0.001, 0.01, 0.1, 1, 10, 100), labels=  c(0.001, 0.01, 0.1, 1, 10, 100))+
  #annotation_logticks(sides="b")+
  #facet_grid(~locality, scales='free')+
  scale_shape_manual(values=c(16, 17, 15, 16, 17), name='')+
  scale_fill_manual(values = c( "#F0E442","#E69F00", "#56B4E9","#e31a1c"), name ="")+
  scale_color_manual(values = c( "#F0E442","#E69F00", "#56B4E9","#e31a1c"), name="") +
  #ggtitle(bquote(alpha ~~ R^2))+
  ylab(expression(R^2))+
  xlab(expression(sigma[epsilon]^2))+
  theme(legend.position=c(1,0), legend.justification=c(1,0), strip.background=element_blank(), legend.background=element_blank())
p3b

corTestPretty(subset(modelFit, kernel=='BMT')$beta ,subset(modelFit, kernel=='BMT')$R2, maxBF=Inf)
corTestPretty(subset(modelFit, kernel=='GP')$beta ,subset(modelFit, kernel=='GP')$R2, maxBF=Inf)
p3c<- ggplot(subset(modelFit, kernel %in% c('GP', 'BMT')), aes(x=beta, y=R2, color=kernel, fill=kernel, shape = kernel))+
  #geom_jitter(alpha=0.4)+
  stat_smooth(method='lm', alpha = 0.1, fullrange=TRUE) +
  geom_point(alpha=0.7)+
  #geom_rug(alpha=0.5)+
  theme_classic()+
  xlim(c(0,max(upperlimit('GP', 'UCB', c('beta')), upperlimit('BMT', 'UCB', c('beta')))))+ #Max axis limit based on the larger of the two Tukey fences
  #scale_x_continuous(trans=log10_trans(), breaks = c(0.001, 0.01, 0.1, 1, 10, 100), labels=  c(0.001, 0.01, 0.1, 1, 10, 100))+
  #annotation_logticks(sides="b")+
  #facet_grid(~locality, scales='free')+
  ylab(expression(R^2))+
  xlab(expression(beta))+
  scale_shape_manual(values=c(16, 17, 15, 16, 17), name='')+
  scale_fill_manual(values = c("#009E73", "#F0E442","#E69F00", "#56B4E9","#e31a1c"), name ="")+
  scale_color_manual(values = c("#009E73", "#F0E442","#E69F00", "#56B4E9","#e31a1c"), name="") +
  theme(legend.position=c(0,0), legend.justification=c(0,0), strip.background=element_blank(), legend.background=element_blank())
p3c

corTestPretty(subset(modelFit, kernel=='GP')$tau ,subset(modelFit, kernel=='GP')$R2, maxBF=Inf)
corTestPretty(subset(modelFit, kernel=='BMT')$tau ,subset(modelFit, kernel=='GP')$R2, maxBF=Inf)
corTestPretty(subset(modelFit, kernel=='dNN')$tau ,subset(modelFit, kernel=='GP')$R2, maxBF=Inf)
corTestPretty(subset(modelFit, kernel=='kNN')$tau ,subset(modelFit, kernel=='GP')$R2, maxBF=Inf)
p3d<- ggplot(subset(modelFit, kernel %in% c('GP', 'BMT', 'dNN', 'kNN')), aes(x=tau, y=R2, color=kernel, fill=kernel, shape = kernel))+
  #geom_jitter(alpha=0.4)+
  stat_smooth(method='lm', alpha = 0.1, fullrange=TRUE) +
  geom_point(alpha=0.7)+
  #geom_rug(alpha=0.5)+
  theme_classic()+
  xlim(c(0,max(upperlimit('GP', 'UCB', c('tau')), upperlimit('BMT', 'UCB', c('tau')), upperlimit('dNN', 'GM', c('tau')), upperlimit('kNN', 'GM', c('tau')))))+ #Max axis limit based on the larger of the two Tukey fences
  #scale_x_continuous(trans=log10_trans(), breaks = c(0.001, 0.01, 0.1, 1, 10, 100), labels=  c(0.001, 0.01, 0.1, 1, 10, 100))+
  #annotation_logticks(sides="b")+
  #facet_grid(~locality, scales='free')+
  ylab('R2')+
  scale_shape_manual(values=c(16, 17, 15, 14), name='')+
  scale_fill_manual(values = c("#009E73", "#F0E442","#E69F00", "#56B4E9","#e31a1c"), name ="")+
  scale_color_manual(values = c("#009E73", "#F0E442","#E69F00", "#56B4E9","#e31a1c"), name="") +
  ylab(expression(R^2))+
  xlab(expression(tau))+
  theme(legend.position=c(1,1), legend.justification=c(1,1), strip.background=element_blank(), legend.background=element_blank())
p3d

corTestPretty(subset(modelFit, kernel=='GP')$omega ,subset(modelFit, kernel=='GP')$R2, maxBF=Inf)
corTestPretty(subset(modelFit, kernel=='BMT')$omega ,subset(modelFit, kernel=='GP')$R2, maxBF=Inf)
corTestPretty(subset(modelFit, kernel=='dNN')$omega ,subset(modelFit, kernel=='GP')$R2, maxBF=Inf)
corTestPretty(subset(modelFit, kernel=='kNN')$omega ,subset(modelFit, kernel=='GP')$R2, maxBF=Inf)
p3e<- ggplot(subset(modelFit, kernel %in% c('GP', 'BMT', 'dNN', 'kNN')), aes(x=omega, y=R2, color=kernel, fill=kernel, shape = kernel))+
  #geom_jitter(alpha=0.4)+
  stat_smooth(method='lm', alpha = 0.1, fullrange=TRUE) +
  geom_point(alpha=0.7)+
  #geom_rug(alpha=0.5)+
  theme_classic()+
  xlim(c(0,max(upperlimit('GP', 'UCB', c('omega')), upperlimit('BMT', 'UCB', c('omega')), upperlimit('dNN', 'GM', c('omega')), upperlimit('kNN', 'GM', c('omega')))))+ #Max axis limit based on the larger of the two Tukey fences
  #scale_x_continuous(trans=log10_trans(), breaks = c(0.001, 0.01, 0.1, 1, 10, 100), labels=  c(0.001, 0.01, 0.1, 1, 10, 100))+
  #annotation_logticks(sides="b")+
  #facet_grid(~locality, scales='free')+
  ylab('R2')+
  scale_shape_manual(values=c(16, 17, 15, 14), name='')+
  scale_fill_manual(values = c("#009E73", "#F0E442","#E69F00", "#56B4E9","#e31a1c"), name ="")+
  scale_color_manual(values = c("#009E73", "#F0E442","#E69F00", "#56B4E9","#e31a1c"), name="") +
  ylab(expression(R^2))+
  xlab(expression(omega))+
  theme(legend.position='right', strip.background=element_blank(), legend.background=element_blank())
p3e


##############################################################################
#Bonus Round plots
##############################################################################
modelDF <- readRDS('modelResults/BanditBonusRoundmodelDF.Rds')
#Mixed effect regression lines (models run in banditBonus.R), but being loaded here
mGPpredictions <- run_model(brm(participantJudgments ~ modelPredictions  + (1+modelPredictions|id) , data = subset(modelDF, model=='GP'), cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'mGPpredictions')
mBMTpredictions <- run_model(brm(participantJudgments ~ modelPredictions  + (1+modelPredictions|id) , data = subset(modelDF, model=='BMT'), cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'mBMTpredictions')
mdNNpredictions <- run_model(brm(participantJudgments ~ modelPredictions  + (1+modelPredictions|id) , data = subset(modelDF, model=='dNN'), cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'mdNNpredictions')
mkNNpredictions <- run_model(brm(participantJudgments ~ modelPredictions  + (1+modelPredictions|id) , data = subset(modelDF, model=='kNN'), cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'mkNNpredictions')

#Predict new data
newdat <- expand.grid(modelPredictions=seq(1,100))
gppreds <- fitted(mGPpredictions, re_formula = NA, newdata = newdat, probs = c(0.025, 0.975))
bmtPreds <- fitted(mBMTpredictions, re_formula = NA, newdata = newdat, probs = c(0.025, 0.975))
dNNPreds <- fitted(mdNNpredictions, re_formula = NA, newdata = newdat, probs = c(0.025, 0.975))
kNNPreds <- fitted(mkNNpredictions, re_formula = NA, newdata = newdat, probs = c(0.025, 0.975))
fixedDF <- data.frame(modelPredictions = rep(newdat$modelPredictions,4),model = rep(c('GP', 'BMT', 'dNN', 'kNN'), each = 100), 
                      participantJudgments = c(gppreds[,1], bmtPreds[,1], dNNPreds[,1], kNNPreds[,1]), 
                      lower = c(gppreds[,3], bmtPreds[,3], dNNPreds[,3], kNNPreds[,3]),  
                      upper = c(gppreds[,4], bmtPreds[,4], dNNPreds[,4], kNNPreds[,4]) )
fixedDF$model <- factor(fixedDF$model, levels=c('GP', 'BMT', 'dNN', 'kNN'))
#Bayesian R2
annText <- data.frame(model = c('GP', 'BMT', 'dNN', 'kNN'), modelPredictions=50, participantJudgments=100, Text = c("paste(italic(R)^2 , \" = .32\")", "paste(italic(R)^2 , \" = .18\")", "paste(italic(R)^2 , \" = .33\")", "paste(italic(R)^2 , \" = .33\")"))
#Add breaks for computing means along 100 equally spaced intervals
modelDF$modelPredictionsBreaks<- cut(modelDF$modelPredictions, breaks=seq(1,100))
modelDF$modelPredictionsBreaks <- as.numeric(modelDF$modelPredictionsBreaks) +.5

pModelPredictions <- ggplot(fixedDF, aes(x=modelPredictions, y=participantJudgments, colour=model, fill = model)) +
  theme_classic()+
  facet_grid(~model)+
  #geom_smooth(data = modelDF, method='lm', fill=NA, aes(group=id),size=.5,color = 'black', alpha = 0.1) +
  stat_summary(data = modelDF, aes(x = modelPredictionsBreaks), fun.y = mean, geom='point',color = 'black', alpha = 0.5) +
  geom_line(data = fixedDF,  size = 1)+ #Fixed efects
  geom_ribbon(data = fixedDF, aes(ymin=lower, ymax = upper), color = NA, alpha = 0.1)+
  #scale_color_manual(values = c( "#009E73", "#F0E442","#e41a1c","#E69F00", "#56B4E9","#e31a1c"), name ="")+
  scale_color_manual(values = c( "#009E73", "#F0E442","#E69F00", "#56B4E9","#e31a1c"), name ="")+
  scale_fill_manual(values = c( "#009E73", "#F0E442","#E69F00", "#56B4E9","#e31a1c"), name ="")+
  coord_cartesian(ylim = c(0,100), xlim = c(0,100))+
  ylab('Participant Judgments')+
  xlab('Model Predictions')+
  geom_text(data = annText, aes(label = Text) , parse=T, color = 'black')+
  theme(legend.position = 'none', legend.key=element_blank(),legend.background=element_blank(), strip.background=element_blank(), text = element_text(size=12,  family="sans"))
pModelPredictions


#Error
errorDF <- ddply(modelDF, c('id', 'model'), summarize, RMSE = sqrt(mean((modelPredictions - participantJudgments)^2)))
randomPreds <- sqrt(mean(sapply(1:10000, FUN=function(i) (runif(n = 1, min=0, max = 100) - sample(bonusDF$judgment, 1))^2)))
modelError <- ggplot(errorDF, aes(x = model, y = RMSE, color = model))+
  geom_hline(yintercept = randomPreds, linetype = 'dashed')+
  geom_quasirandom(alpha = 0.5)+
  stat_summary(fun.data = mean_cl_boot, geom='errorbar', color= 'black', width = 0.1)+
  stat_summary(fun.y = mean, geom='point', fill = 'white', color= 'black', shape = 23)+
  scale_color_manual(values = c( "#009E73", "#F0E442","#E69F00", "#56B4E9","#e31a1c"), name ="")+
  theme_classic()+
  coord_cartesian(ylim=c(5, 65))+
  scale_y_continuous(breaks = pretty_breaks(n=3)(5:60))+
  xlab('')+
  geom_signif(comparisons=list( c("GP", "BMT"), c("GP", "dNN"), c("GP", "kNN"), c("dNN", "kNN") ), annotations=c("italic(BF)[10] == 3.5 %*% 10^6","italic(BF)[10] == 0.11", "italic(BF)[10] == 0.13",  "italic(BF)[10] == 0.14"), parse=T, y_position = c(40,50,60, 40), tip_length = 0, col="black")+
  theme(legend.position='none', strip.background=element_blank(), legend.key=element_blank(),legend.background = element_rect(fill=NA),text = element_text(size=12,  family="sans"))
modelError

#T-tests
ttestPretty(subset(errorDF, model == 'GP')$RMSE, subset(errorDF, model == 'BMT')$RMSE, paired=T, maxBF=Inf)
ttestPretty(subset(errorDF, model == 'GP')$RMSE, subset(errorDF, model == 'dNN')$RMSE, paired=T)
ttestPretty(subset(errorDF, model == 'GP')$RMSE, subset(errorDF, model == 'kNN')$RMSE, paired=T)
ttestPretty(subset(errorDF, model == 'dNN')$RMSE, subset(errorDF, model == 'kNN')$RMSE, paired=T)

#Participants best predicted
library(nnet)
bestFit <- table(sapply(unique(errorDF$id), FUN=function(i){
  subd <- subset(errorDF, id==i)
  subd$model[which.is.max(subd$RMSE*-1)] #convert RMSe to negative to use which.is.max
}))

bestFitDF <- data.frame(model = c('GP', 'BMT', 'dNN', 'kNN'), bestFit =bestFit) #convert to dataframe
bestFitDF$model <- factor(bestFitDF$model, levels = c('GP', 'BMT', 'dNN', 'kNN'))

pBestFit <- ggplot(bestFitDF, aes(x = model, y = bestFit, fill = model))+
  geom_bar(stat='identity', color = 'black')+
  theme_classic()+
  xlab('')+
  ylab('Participants Best Described\n (Bonus Round)')+
  scale_fill_manual(values = c( "#009E73", "#F0E442","#E69F00", "#56B4E9","#e31a1c"), name ="")+
  theme(legend.position='none')
pBestFit  


#Mixed effects model
gpConf <- run_model(brm(participantConfidence ~ modeluncertainty + (1+modeluncertainty|id) , data = subset(modelDF, model=='GP'), cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'mGPConf')
gpConfNull <- run_model(brm(participantConfidence ~ 1 + (1+modeluncertainty|id) , data = subset(modelDF, model=='GP'), cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'gpConfNull')
#bayes_R2(gpConf)# 0.4980354 
bayes_factor(gpConf, gpConfNull) #453260.16300
#Now try to predict 
fixedTerms <- fixef(gpConf)

#Now generate predictions, removing id as a random effect
xseq <- seq(0,1, length.out = 100)
newdat <-data.frame(modeluncertainty = xseq)
preds <- fitted(gpConf, re_formula = NA, newdata = newdat, probs = c(0.025, 0.975))
#create new fixed effects dataframe
fixedDF <- data.frame( modeluncertainty =xseq, participantConfidence = preds[,1], lower = preds[,3], upper = preds[,4] )

modelDF$uncertaintyBreaks<- cut(modelDF$modeluncertainty, breaks= seq(0,1, length.out=100))
levels(modelDF$uncertaintyBreaks) <- seq(0,1, length.out=100)
modelDF$uncertaintyBreaks <- as.numeric(as.character(modelDF$uncertaintyBreaks))

modelconf <- ggplot(subset(modelDF, model %in% c('GP', 'BMT')), aes(x = participantConfidenceRank, y = modelUncertaintyRank, color = model))+
  #geom_point(alpha=0.2)+
  stat_summary(fun.y = mean, geom= 'point', size = 2)+
  stat_summary(data = subset(modelDF, model == 'GP'), aes(group=model),fun.data = mean_cl_boot, geom= 'errorbar', color = 'black')+
  geom_smooth(se=F, method = 'lm', alpha = 0.8)+
  scale_color_manual(values =c( "#009E73", "#F0E442"), name ='')+
  xlab('Participant Confidence Rank')+
  ylab('Model Uncertainty Rank')+
  theme_classic()+
  scale_x_continuous(breaks=c(0,2,4, 6,8,10))+
  scale_y_continuous(breaks=c(0,2,4, 6,8,10))+
  coord_cartesian(ylim = c(1,10), xlim = c(1,10), expand = c(0.5))+
  annotate("text", x = 5, y = 10, label = "paste(italic(b)[gpUncertainty] , \" = -3.1, 95% HPD: [-4.4, -1.9]\")", parse = TRUE)+
  #annotate("text", x = 5, y = 10, label = "paste(italic(b)[gpUncertainty] , \" = -3.1\", \", \",  italic(BF)[10] == 4.5 %*% 10^5)", parse = TRUE)+
  theme(legend.position='None', strip.background=element_blank(), legend.key=element_blank(),legend.background = element_rect(fill=NA),text = element_text(size=12,  family="sans"))
modelconf

#TAble
tab_model(mGPpredictions,mBMTpredictions, mdNNpredictions, mkNNpredictions, gpConf)

########################################################################
# Put it all together
########################################################################
#Main Plot
topRow <- cowplot::plot_grid(trialPlot,pRepeatReward, pDistanceReward, rows = 1, labels = 'auto') #from banditExpBehaviorPlots.R
middleRow <- cowplot::plot_grid(p1b, pLearningCurve, pBestFit, rel_widths = c(1,1,1),nrow = 1,  labels = c('d','e','f'))
bottom <- cowplot::plot_grid(pModelPredictions, modelconf, rel_widths = c(2,1), labels = c('g','h'))
mainPlot <-  cowplot::plot_grid(topRow, middleRow, bottom, ncol=1, labels = c('','')) #TODO add bonus round
mainPlot
ggsave('plots/exp2Results.pdf', mainPlot, width = 12, height = 9, units = 'in')

#Save parameters
ggsave('plots/exp2Params.pdf', pParams, width = 10, height = 6, units = 'in')

#SI Plots
modelDiffPlots <- cowplot::plot_grid(pModelDiffGP_BMT, pModelDiffGP_dNN, pModelDiffGP_kNN,pModelDiffGP_Sticky, nrow = 2, labels = 'auto')
modelDiffPlots
ggsave('plots/modelDiffExp2.pdf', modelDiffPlots, width = 10, height = 8, units = 'in')


pParamR2Top <- cowplot::plot_grid(p3a, p3b, p3c, nrow = 1, labels='auto')
pParamR2Bottom <- cowplot::plot_grid(p3d, p3e, rel_widths = c(1.3, 1.7), labels=c('d','e'))
pParamR2 <- cowplot::plot_grid(pParamR2Top, pParamR2Bottom, ncol = 1)
ggsave('plots/ParamR2.pdf', pParamR2, width = 12, height = 8, units = 'in')

