#BanditBonus.R
rm(list=ls())

packages <- c('tidyr',"jsonlite", "ggplot2", "scales", 'plyr',  "brms", "RColorBrewer", 'lmerTest', 'ggridges', 'cowplot', 'brms', "igraph", 'matrixcalc', 'Matrix', 'dplyr', 'corrplot', 'pander', 'BayesFactor', 'ggbeeswarm', 'viridis', 'lsr', 'ggsignif')
lapply(packages, require, character.only = TRUE)

source('exportImportGraph.R')
source('models.R')
source('utilities.R')
source('statisticalTests.R')
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colfunc<-colorRampPalette(c("#0D0887FF", "#CC4678FF", "#F0F921FF"))

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
eigenCentrality = list()

for (i in 1:length(edges)){
  g <- graph_from_edgelist(as.matrix(edges[[i]]), directed=F) 
  V(g)$reward  <- nodes[[i]]$reward
  V(g)$x <- nodes[[i]]$x / 5 #divide by scaling factor used in graph generation for visualization purposes
  V(g)$y <- nodes[[i]]$y / 5 
  #plot(g, vertex.size = V(g)$reward)
  graphList[[i]] <-g
  distanceTables[[i]] <- distances(g)
  nodeDegrees[[i]] <- igraph::degree(g)
  eigenCentrality[[i]] <- eigen_centrality(g)$vector
}


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
####################################################################################
# Participant data
####################################################################################
#Set of data that needs to be recovered from experiment

df <- dataImportBandit('../data/banditTask/networkBandit.csv', delim=',')
df$id <- factor(df$id)


####################################################################################
# Bonus data
####################################################################################

rawdata <- read.csv('../data/banditTask/networkBandit.csv', sep=',', header = T)
#remove empty rows
rawdata <- rawdata[!grepl("NULL",rawdata$task_end) & !grepl("NULL",rawdata$task_end),] 
#extract sampleSize, records, and dat (deeply nested json data)
sampleSize <-1:nrow(rawdata)
cleanedList <- sampleSize[!sampleSize %in% c(45, 92)] #Exclude participants with Missing data
bonusDF <- data.frame()
#Loop through data and build df
for (i in cleanedList){
  dat = rawdata[i,]
  id <- dat$id
  experimentData <- fromJSON(as.character(dat$experimentData), simplifyVector=T, simplifyMatrix=T)
  g <- graphList[[experimentData$envOrder[15]+ 1]] #plus one is to switch from base_0 to base_1 indexing
  nodeDistances <- distances(g)
  bonusDat <- experimentData$bonusCollect$bonusStimuli
  chosen <- experimentData$xcollect[10,22]
  choice <- rep(F, 10)
  choice[which(bonusDat$x==chosen)] <- T
  RTs <- c(NA,bonusDat$ts[2:10] -  bonusDat$ts[1:9])
  #some graph statistics
  nodeDegree <- degree(g)[bonusDat$x]
  observedNodes <- experimentData$xcollect[10,1:21]
  distanceOtherNodes <- nodeDistances[,bonusDat$x][observedNodes,] #Calculate distances of bonus nodes to observed nodes
  closestObservation <- sapply(1:10, FUN=function(n) min(distanceOtherNodes[n,]))
  trueValue <-normalizeMinMax(V(g)$reward[bonusDat$x], experimentData$minList[15], experimentData$maxList[15]) #True value, rescaled based on the min and max rescaling that participants also observed
  #build dataframe
  dummy <-data.frame(id = rep(id, 10), nodeId = bonusDat$x, nodeDegree= nodeDegree,closestObservation=closestObservation, judgment = bonusDat$givenValue, trueValue = trueValue, confidence = bonusDat$howCertain, RT = RTs, chosen = choice , graphId =experimentData$envOrder[15] + 1 )
  bonusDF<-rbind(bonusDF, dummy)
}

bonusDF$RSE <- sqrt((bonusDF$trueValue - bonusDF$judgment)^2)
#Compare to random performance
randomGuesses <- runif(10000,0,100)
trueValues <- sample(bonusDF$trueValue, 10000,replace = T) #sample with replacement
error <- sapply(1:10000, FUN=function(i) randomGuesses[i] - trueValues)
randomPerformance <- sqrt(mean(error^2))

#Add in the eigen centrality of each node
eigenCentralityVec <- rep(NA, nrow(bonusDF))
for (i in 1:nrow(bonusDF)){
  subd <- bonusDF[i,]  
  eigenCentralityVec[i] <- eigenCentrality[[subd$graphId]][subd$nodeId]
}
bonusDF$eigenCentrality <- eigenCentralityVec

####################################################################################
# Plot behavioral results
####################################################################################
#average by participant
pDF <- ddply(bonusDF, ~id, plyr::summarize, RMSE = sqrt(mean((trueValue -judgment)^2)), meanConfidence = mean(confidence))
pDF$score <- ddply(subset(df, round<10), ~id,  plyr::summarize, performance = mean(rewardObjective))$performance

#Bandit performance to judgment performance
corTestPretty(pDF$score, pDF$RMSE)

pPerJudg <- ggplot(as.data.frame(pDF), aes(x=score, y = RMSE))+
  geom_hline(yintercept = randomPerformance, linetype = 'dashed')+
  geom_point(alpha = 0.8, color = '#377eb8')+
  theme_classic()+
  xlab('Mean Bandit Score')+
  geom_smooth(method='lm', color = "#e41a1c", fill = "#e41a1c")+
  #geom_smooth(method = 'lm', se = FALSE, color= "#e41a1c")+
  annotate("text", x = 60, y = 50, label = "paste(italic(r) , \" = -.14\", \", \", italic(BF[10]) ,\" = 0.55\" )", parse = TRUE)
pPerJudg


#judgment performance to accuracy
corTestPretty(bonusDF$confidence, bonusDF$RSE) 
#Mixed models
bonusJudgmentRMSE <- run_model(brm(RSE ~ confidence + (1 + confidence|id), data =bonusDF,  cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'bonusJudgmentRMSE')
bonusJudgmentRMSENull <- run_model(brm(RSE ~ 1 + (1 + confidence|id), data =bonusDF,  cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'bonusJudgmentRMSENull')
#bayes_factor(bonusJudgmentRMSE, bonusJudgmentRMSENull) #0.76379
#bayes_R2(bonusJudgmentRMSE) #slow
fixedTerms <- fixef(bonusJudgmentRMSE)#Look at fixed terms

#Now generate predictions, removing id as a random effect
xseq <- seq(1,11)
newdat <-data.frame(confidence = xseq)
preds <- fitted(bonusJudgmentRMSE, re_formula = NA, newdata = newdat, probs = c(0.025, 0.975))
#create new fixed effects dataframe
fixedDF <- data.frame( confidence =xseq, RSE = preds[,1], lower = preds[,3], upper = preds[,4] )

pbonusJudgmentRMSE<- ggplot(as.data.frame(bonusDF), aes(x=confidence, y = RSE))+
  #geom_hline(yintercept = randomPerformance, linetype = 'dashed')+
  #geom_quasirandom(alpha = 0.8, color = '#377eb8')+
  #stat_smooth(aes(group=id), method = 'lm', geom='line', alpha=0.5, se=FALSE,color='#377eb8')+
  stat_summary(fun.y = mean, geom='point', color ='#377eb8')+
  stat_summary(fun.data = mean_cl_boot, geom='errorbar', color ='#377eb8', width = 0)+
  geom_line(data = fixedDF,  size = 1, color = "#e41a1c")+ #Fixed efects
  geom_ribbon(data = fixedDF, aes(ymin=lower, ymax = upper), color = NA, alpha = 0.4, fill = "#e41a1c" )+
  theme_classic()+
  labs(x='Confidence', y = 'Root Squared Error')+
  annotate("text", x = 6, y = 35, label = "paste(italic(b)[confidence] , \" = -0.1\", \", \", italic(BF)[10] ,\" = 0.8\" )", parse = TRUE)
pbonusJudgmentRMSE

corTestPretty(pDF$meanConfidence, pDF$RMSE) #average per individual
p2agg<- ggplot(as.data.frame(pDF), aes(x=meanConfidence, y = RMSE))+
  geom_hline(yintercept = randomPerformance, linetype = 'dashed')+
  geom_quasirandom(alpha = 0.8, color = '#377eb8')+
  theme_classic()+
  xlab('Mean Confidence')+
  geom_smooth(method = 'lm', color= "#e41a1c", fill= "#e41a1c")+
  annotate("text", x = 6, y = 50, label = "paste(italic(r) , \" = .19\", \", \", italic(BF)[10] ,\" = 1.3\" )", parse = TRUE)
p2agg

#Judgment accuracy to degree
corTestPretty(bonusDF$nodeDegree, bonusDF$RSE) #average per individual
p3a<- ggplot(as.data.frame(bonusDF), aes(x=nodeDegree, y = RSE))+
  geom_hline(yintercept = randomPerformance, linetype = 'dashed')+
  geom_quasirandom(alpha = 0.8, color = '#377eb8')+
  theme_classic()+
  geom_smooth(method = 'lm', se = FALSE, color= "#e41a1c")+
  annotate("text", x = 3, y = 85, label = "paste(italic(r) , \" = -.09\", \", \", italic(BF) ,\" = 4.7\" )", parse = TRUE)
p3a

#Judgment accuracy to eigenCentrality
corTestPretty(bonusDF$eigenCentrality, bonusDF$RSE) #average per individual
p3ab<- ggplot(as.data.frame(bonusDF), aes(x=eigenCentrality, y = RSE))+
  geom_hline(yintercept = randomPerformance, linetype = 'dashed')+
  geom_quasirandom(alpha = 0.2, color = '#377eb8')+
  theme_classic()+
  geom_smooth(method = 'lm', se = FALSE, color= "#e41a1c")+
  annotate("text", x = .5, y = 80, label = "paste(italic(r) , \" = .05\", \", \", italic(BF)[10] ,\" = 0.32\" )", parse = TRUE)
p3ab


corTestPretty(bonusDF$nodeDegree, bonusDF$confidence) #average per individual
p3b<- ggplot(as.data.frame(bonusDF), aes(x=nodeDegree, y = confidence))+
  #geom_hline(yintercept = randomPerformance, linetype = 'dashed')+
  geom_quasirandom(alpha = 0.8, color = '#377eb8')+
  geom_smooth(method = 'lm', se = FALSE, color= "#e41a1c")
  #annotate("text", x = .5, y = 75, label = "paste(italic(r) , \" = .05\", \", \", italic(BF)[10] ,\" = 0.32\" )", parse = TRUE)
p3b


#Judgment confidence to RT
corTestPretty(bonusDF$confidence, log(bonusDF$RT)) #average per individual
p4<- ggplot(as.data.frame(bonusDF), aes(x=RSE, y =log10(RT)))+
  #geom_hline(yintercept = randomPerformance, linetype = 'dashed')+
  geom_quasirandom(alpha = 0.8, color = '#377eb8')+
  geom_smooth(method = 'lm', se = FALSE, color= "#e41a1c")
  #annotate("text", x = 3, y = 85, label = "paste(italic(r) , \" = -.09\", \", \", italic(BF) ,\" = 5\" )", parse = TRUE)
p4

##################################################################################
####Main Plots
################################################################################

#How well did participants perform?
ttestPretty(pDF$RMSE, mu=randomPerformance, maxBF=Inf) #Error was better than random
corTestPretty(bonusDF$judgment, bonusDF$trueValue, maxBF=Inf) #Judgments were correlated

corTestPretty(subset(bonusDF, chosen==T)$judgment, subset(bonusDF, chosen==T)$trueValue) #Chosen judgments
corTestPretty(subset(bonusDF, chosen==F)$judgment, subset(bonusDF, chosen==F)$trueValue) #Chosen judgments

#Compare RMSE of chosen judgment vs. non-chosen
chosenRSE <- subset(bonusDF, chosen==T)$RSE
nonChosenRMSE <- ddply(subset(bonusDF, chosen==F), ~id, plyr::summarize, RMSE = mean(RSE))
ttestPretty(chosenRSE, nonChosenRMSE$RMSE, paired=T) #no difference between accuracy of chosen and non-chosen

bonusDF$nodeType <- "Not Chosen"
bonusDF[bonusDF$chosen==T, 'nodeType'] <- "Chosen"
bonusDF$nodeType <- factor(bonusDF$nodeType)

#Correlation
pJudgment <- ggplot(as.data.frame(bonusDF), aes(x=judgment, y = trueValue))+
  geom_point(alpha = 0.8, color = '#377eb8')+
  geom_smooth(method = 'lm', color= "#e41a1c")+
  coord_cartesian(ylim = c(0,100), xlim = c(0,100), expand = c(0.5))+
  facet_grid(~nodeType)+
  theme_classic()+
  xlab('Bonus Judgment')+
  ylab('True Value')+
  theme(legend.position='none', strip.background=element_blank(), legend.key=element_blank(),legend.background = element_rect(fill=NA))
pJudgment

#Judgment Error of chosen vs. not chosen
ErrorChosenDF <- ddply(bonusDF, ~id+chosen, plyr::summarize, RMSE = mean(RSE))
ErrorChosenDF$nodeType <- ifelse(ErrorChosenDF$chosen, 'Chosen', 'Not Chosen')
ttestPretty(subset(ErrorChosenDF, chosen==T)$RMSE, subset(ErrorChosenDF, chosen==F)$RMSE, paired=T)

pErrorChosenDF <- ggplot(ErrorChosenDF, aes(x=nodeType, y = RMSE))+
  geom_hline(yintercept = randomPerformance, linetype = 'dashed')+
  geom_quasirandom(alpha=0.8, color="#377eb8")+
  stat_summary(fun.data = mean_cl_boot, geom='errorbar', width = 0, color = "#e41a1c" )+
  stat_summary(fun.y = mean, geom='point', color = '#e41a1c', size = 2)+
  xlab('')+
  ylab('Judgment Error (RMSE)')+
  theme_classic()+
  ylim(c(0,75))+
  #scale_y_continuous(breaks=seq(1, 11, length.out = 3))+
  #scale_color_manual(values = c("#d95f02", '#1b9e77' ))+
  theme(legend.position='none', strip.background=element_blank(), legend.key=element_blank(),legend.background = element_rect(fill=NA))+
  geom_signif(comparisons = list(c('Chosen', 'Not Chosen')), annotation = c("italic(BF)[10] == 0.59"),parse=T, tip_length = 0)
pErrorChosenDF

#Confidence and chosen vs. not chosen
confDF <- ddply(bonusDF, ~id+chosen, plyr::summarize, meanConfidence = mean(confidence))
ttestPretty(subset(confDF, chosen==T)$meanConfidence, subset(confDF, chosen==F)$meanConfidence, paired=T)
confDF$nodeType <- ifelse(confDF$chosen, 'Chosen', 'Not Chosen')
pConfidence <- ggplot(confDF, aes(x=nodeType, y = meanConfidence))+
  geom_quasirandom(alpha=0.8, color="#377eb8")+
  stat_summary(fun.data = mean_cl_boot, geom='errorbar', width = 0, color = "#e41a1c" )+
  stat_summary(fun.y = mean, geom='point', color = '#e41a1c', size = 2)+
  xlab('')+
  ylab('(Mean) Confidence')+
  theme_classic()+
  scale_y_continuous(breaks=seq(1, 11, length.out = 3), lim=c(1,12))+
  #scale_color_manual(values = c("#d95f02", '#1b9e77' ))+
  theme(legend.position='none', strip.background=element_blank(), legend.key=element_blank(),legend.background = element_rect(fill=NA))+
  geom_signif(comparisons = list(c('Chosen', 'Not Chosen')), annotation=c("italic(BF)[10] == 21"),parse=T, tip_length = 0)
pConfidence 

####################################################################################
# Regression for choice
####################################################################################

#Bayes factor comparison
mChoice  <- run_model(brm(chosen ~ judgment +confidence + ( 1+ judgment +confidence| id ), data =bonusDF, family = 'bernoulli', cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'mChoice')
mChoiceNoJudgment  <- run_model(brm(chosen ~ 1 + confidence + ( 1+ judgment +confidence| id ), data =bonusDF, family = 'bernoulli', cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'mChoiceNoJudgment')
mChoiceNoConf  <- run_model(brm(chosen ~ confidence + 1 + ( 1+ judgment +confidence| id ), data =bonusDF, family = 'bernoulli', cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'mChoiceNoConf')
#bayes_factor(mChoice, mChoiceNoJudgment) #Effect of judgment #BF=9456008083261.51172
#bayes_factor(mChoice, mChoiceNoConf) #Effect of confidence #BF= 10240207482649.02148
#bayes_R2(mChoice) #0.1037533
fixedTerms <- fixef(mChoice)#Look at fixed terms
summary(mChoice)

library(sjPlot)
pMChoice <- plot_model(mChoice,vline.color = NA, bpe.color = "#e41a1c", bpe.style='point')+ theme_classic() + geom_hline(yintercept = 1, linetype = 'dashed')+
  ylim(c(0.95,1.2)) + scale_fill_manual(values = c('#377eb8')) + scale_color_manual(values = c('#377eb8')) + ggtitle('Bonus Round Choice')
pMChoice

library(tidybayes)
logOdds <- mChoice %>% gather_draws(b_judgment,b_confidence ) %>% median_hdi()
print(exp(logOdds[c(".value",".lower",".upper")])) #Odds ratio
####################################################################################
# Model based results
####################################################################################
source('models.R')
library('DescTools')

#parameter estimates
bmtPars <-  read.csv('rationalModels/parameters/BMT.csv')
gpPars <- read.csv('rationalModels/parameters/DF.csv')
knnPars <-  read.csv('rationalModels/parameters/kNN.csv')
dnnPars <-  read.csv('rationalModels/parameters/dNN.csv')

modelDF <- data.frame()

for (pid in unique(df$id)){
  subdf <- subset(bonusDF, id==pid)
  subd <- subset(df, id==pid)
  g <- graphList[[subdf$graphId[1]]]
  distanceMatrix <- distanceTables[[subdf$graphId[1]]]
  #compute model predictions
  observedNodes <- subset(subd, round == 10)$nodeId[1:21]
  observedRewards <- round(subset(subd, round == 10)$rewardSubjective[1:21])
  #BMT always predicts the same value
  bmtpreds <- rep(50,10)
  bmtConf <- rep(5,10)
  #GP-DF predictions
  K<-diffusionKernel(g,alpha= median(subset(gpPars, participant == pid)$alpha))
  gpPost <- gpr(X.test = 1:64, X = observedNodes, Y = observedRewards, mu_0 = 50, k=K)
  gppreds <- floorCeiling(gpPost[subdf$nodeId,'mu'], min=0, max=100)
  gpConf <- gpPost[subdf$nodeId,'var']
  #kNN preds
  kNNPreds <- floorCeiling(knnChoice(distanceMatrix, X = observedNodes, Y = observedRewards ,mu_0 = 50, kValue  = median(subset(knnPars, participant == pid)$k))$mu[subdf$nodeId], min=0, max = 100)
  #dNN preds
  dNNPreds <- floorCeiling(dnnChoice(distanceMatrix, X = observedNodes, Y = observedRewards , mu_0 = 50,  dist  = median(subset(dnnPars, participant == pid)$d))$mu[subdf$nodeId], min=0, max = 100)
  #put it in dataframe
  dummy <- data.frame(id = rep(pid, 40), modelPredictions = c(bmtpreds, gppreds, kNNPreds, dNNPreds), modeluncertainty = c(bmtConf, gpConf, rep(NA, 20)),modelUncertaintyRank = c(rep(5.5,10), rank(gpConf), rep(NA, 20)), model = rep(c("BMT", "DF",  'kNN', 'dNN'), each=10), 
                      participantJudgments = rep(subdf$judgment,4), participantConfidence = rep(subdf$confidence,4), participantConfidenceRank = rep(rank(subdf$confidence),4 ))
  modelDF <- rbind(modelDF, dummy)
}

#Judgments
corTestPretty(subset(modelDF, model=='DF')$modelPredictions,subset(modelDF, model=='DF')$participantJudgments)
corTestPretty(subset(modelDF, model=='dNN')$modelPredictions,subset(modelDF, model=='dNN')$participantJudgments)
#correlation test
corDF <- ddply(modelDF, c('id', 'model'), summarize, modelCor = cor(modelPredictions, participantJudgments))
corDF[is.na(corDF$modelCor), 'modelCor']<- 0
corDF$Z <- FisherZ(corDF$modelCor)
ddply(corDF, c('model'), summarize, cor = mean(modelCor, na.rm=T))
ttestPretty(subset(corDF, model=='DF')$Z, subset(corDF, model=='dNN')$Z, paired=T)
ttestPretty(subset(corDF, model=='DF')$Z, subset(corDF, model=='kNN')$Z, paired=T)

#######Reorganize########
modelDF <- subset(modelDF, model %in% c('BMT', 'DF', 'dNN', 'kNN')) #Filter
modelDF$model <- factor(modelDF$model, levels = c('DF', 'BMT', 'dNN', 'kNN')) #reorganize
levels(modelDF$model)[1]<- 'GP' #rename

saveRDS(modelDF, 'modelResults/BanditBonusRoundmodelDF.Rds')
#########################################################################################################################################
#Plot Model predictions against participant judgments
#########################################################################################################################################
#Mixed effect regression lines
#How well are participant judgments explained by the GP predictions?
mGPpredictions <- run_model(brm(participantJudgments ~ modelPredictions  + (1+modelPredictions|id) , data = subset(modelDF, model=='GP'), cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'mGPpredictions')
#bayes_R2(mGPpredictions) # 0.3152814
#gpLoo <- loo(mGPpredictions)
mBMTpredictions <- run_model(brm(participantJudgments ~ modelPredictions  + (1+modelPredictions|id) , data = subset(modelDF, model=='BMT'), cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'mBMTpredictions')
#bayes_R2(mBMTpredictions) # 0.1847933
#bmtLoo <- loo(mBMTpredictions)
mdNNpredictions <- run_model(brm(participantJudgments ~ modelPredictions  + (1+modelPredictions|id) , data = subset(modelDF, model=='dNN'), cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'mdNNpredictions')
#bayes_R2(mdNNpredictions) #0.3303634
#dNNLoo <- loo(mdNNpredictions)
mkNNpredictions <- run_model(brm(participantJudgments ~ modelPredictions  + (1+modelPredictions|id) , data = subset(modelDF, model=='kNN'), cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'mkNNpredictions')
#bayes_R2(mdNNpredictions) # 0.3303634 
#kNNLoo <- loo(mkNNpredictions)
#Compare the different models
#loo_compare(gpLoo,dNNLoo)
#loo_compare(gpLoo,kNNLoo)
#loo_model_weights(list(gpLoo,bmtLoo, dNNLoo, kNNLoo)) #0.218 , 0.000 ,  0.369 , 0.414
#tab_model(mGPpredictions, mBMTpredictions, mdNNpredictions, mkNNpredictions)

#Can we improve the GP model by adding dNN predictions?
wideModelDF <- subset(modelDF, model=='GP')
wideModelDF<- rename(wideModelDF, 'gpPredictions' = 'modelPredictions')#change to less confusing naming
wideModelDF$dNNPredictions <- subset(modelDF, model=='dNN')$modelPredictions #add in dNN predictions
mGP_dNNpredictions <- run_model(brm(participantJudgments ~ gpPredictions + dNNPredictions + (1+gpPredictions+dNNPredictions|id) , data = wideModelDF, cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'mGP_dNNpredictions')
#bayes_R2(mGP_dNNpredictions) #0.3445317
#combinedLoo <- loo(mGP_dNNpredictions)
#loo_compare(combinedLoo,gpLoo)
#bayes_factor(mGP_dNNpredictions, mGPpredictions)
fixef(mGP_dNNpredictions)
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


#Main plot
pModelPredictions <- ggplot(fixedDF, aes(x=modelPredictions, y=participantJudgments, colour=model, fill = model)) +
  theme_classic()+
  facet_grid(~model)+
  #geom_point(data = modelDF, size=.5,color = 'black', alpha = 0.1) +
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
ttestPretty(subset(errorDF, model == 'BMT')$RMSE, mu=randomPreds, maxBF=Inf)
ttestPretty(subset(errorDF, model == 'GP')$RMSE, subset(errorDF, model == 'dNN')$RMSE, paired=T)
ttestPretty(subset(errorDF, model == 'GP')$RMSE, subset(errorDF, model == 'kNN')$RMSE, paired=T)
ttestPretty(subset(errorDF, model == 'dNN')$RMSE, subset(errorDF, model == 'kNN')$RMSE, paired=T)


#######
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
  ylab('Participants Best Described')+
  scale_fill_manual(values = c( "#009E73", "#F0E442","#E69F00", "#56B4E9","#e31a1c"), name ="")+
  theme(legend.position='none')
pBestFit  

##### Confidence Judgments #####

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
  theme(legend.position='None', strip.background=element_blank(), legend.key=element_blank(),legend.background = element_rect(fill=NA),text = element_text(size=12,  family="sans"))
modelconf

  

########################################################################
#Put it together
########################################################################
library(cowplot)

#top <- cowplot::plotgrid(pBestFit,modelconf, pModelPredictions, rel_widths = c(1,4), labels = c('d','e'))

#siPlot <- cowplot::plot_grid(pbonusJudgmentRMSE, pPerJudg, p2agg, pErrorChosenDF, pConfidence, pMChoice, labels = 'auto')
siPlot <- cowplot::plot_grid(pErrorChosenDF, pConfidence, pMChoice, nrow = 1,  labels = 'auto')
siPlot
ggsave('plots/banditBonusRound.pdf', siPlot, width = 12, height = 3, units = 'in')

