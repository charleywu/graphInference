#Exp1 Modeling results
rm(list=ls())

packages <- c('tidyr','dplyr', 'corrplot', 'cowplot', 'sjPlot', 'lmerTest','brms', "jsonlite", 'fields', 'MASS', 'GGally', "ggplot2", "RColorBrewer", "igraph", 'matrixcalc', 'Matrix', 'dplyr', 'corrplot', 'pander', 'BayesFactor', 'ggbeeswarm', 'viridis', 'lsr', 'ggsignif', 'cowplot')
lapply(packages, require, character.only = TRUE)

source('exportImportGraph.R')
source('models.R')
source('utilities.R')
source('statisticalTests.R')




####################################################################################
# Load pre-generated graphs
####################################################################################
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
# BMRS wrapper
####################################################################################

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
# Log likelihood function
####################################################################################
MSEtoLogLikelihood <- function(predictions, trueValues, sigma = 1){
  if (is.null(sigma)){
    sigma = sd(predictions - trueValues) #estimate sigma as the standard deviation of model errors
  }
  Z <- 1 / (sigma * sqrt(2*pi))
  nLL = -(sum(-(predictions -trueValues)^2 / (2 * sigma^2)) - log(Z))
  return(nLL)
}

#predictions <- subset(df, id == 66)$gpDFpred
#trueValues <- subset(df, id == 66)$predictions
###############################################################################################################################################
#Load data
###############################################################################################################################################

#df <- readRDS('data/Exp1/MLE.RDS') #MLE results from `Exp1MLE.R`
#df <- readRDSfolder('modelResults/Exp1/cogsci/') #Cross validated results from `Exp1ModelCV.R``

df <- readRDSfolder('modelResults/Exp1/fullPaper/') #Cross validated results Exp1ModelCV.R`` (includes additional graph operators)
modelNames <- c('GP-DF', 'dNN', 'kNN')

#Process data
participantdf <- df %>% group_by(id) %>% 
  summarize(RMSE_kNN = sqrt(mean((predictions - kNNpred)^2)), MAE_kNN = mean(abs(predictions - kNNpred)), SSE_kNN =  sum((predictions - kNNpred)^2), nLL_kNN= MSEtoLogLikelihood(kNNpred, predictions),
  RMSE_dNN = sqrt(mean((predictions -dNNpred)^2 )), MAE_dNN = mean(abs(predictions - dNNpred)), SSE_dNN =  sum((predictions - dNNpred)^2), nLL_dNN = MSEtoLogLikelihood(dNNpred, predictions),
  RMSE_gpDF = sqrt(mean((predictions -gpDFpred)^2 )),MAE_gpDF = mean(abs(predictions - gpDFpred)), SSE_gpDF =  sum((predictions - gpDFpred)^2), nLL_gpDF= MSEtoLogLikelihood(gpDFpred, predictions),
  alpha= mean(alpha), kNNkValue = mean(k), dNNdelta = mean(delta), participantError = sqrt(mean(errorArray^2)))


#prepare plotting df
modelplottingDF <- data.frame(id = rep(participantdf$id, length(modelNames)), modelName = rep(modelNames, each=nrow(participantdf)), 
                              RMSE = c(participantdf$RMSE_gpDF,participantdf$RMSE_dNN, participantdf$RMSE_kNN), 
                              MAE = c(participantdf$MAE_gpDF, participantdf$MAE_dNN, participantdf$MAE_kNN),
                              nLL = c(participantdf$nLL_gpDF, participantdf$nLL_dNN, participantdf$nLL_kNN),
                              SSE = c(participantdf$SSE_gpDF, participantdf$SSE_dNN, participantdf$SSE_kNN))
modelplottingDF$modelName <- factor(modelplottingDF$modelName, levels=modelNames)

#FIlTER MODELS AND REFACTOR
modelplottingDF <- subset(modelplottingDF, modelName %in% c('GP-DF', 'dNN', 'kNN')) #Omit some models
modelplottingDF$modelName <- factor(modelplottingDF$modelName) #refactor
levels(modelplottingDF$modelName) <- c('GP', 'dNN', 'kNN') # rename

###############################################################################################################################################
#PLOTS
###############################################################################################################################################

#number of participants best fit
bestModel <- modelplottingDF %>% group_by(id) %>% count(modelName, RMSE) %>% filter (RMSE==min(RMSE))
table(bestModel$modelName) 
modelplottingDF$bestModel <- FALSE #set default to zero
for (pid in bestModel$id){
  modelplottingDF[modelplottingDF$id == pid & modelplottingDF$modelName == bestModel[bestModel$id==pid, 'modelName'][[1]], 'bestModel'] <- TRUE
}


pBestFit <- ggplot(subset(modelplottingDF, bestModel==T), aes(x =modelName, fill = modelName, color = modelName))+
  geom_bar(color = 'black')+
  ylab('Participants Best Fit')+
  xlab('')+
  theme_classic()+
  scale_fill_manual(values =c( "#009E73", "#E69F00", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), name = '')+
  theme(legend.position = 'none')
pBestFit

#NLL
randomNLL <-mean(sapply(1:1000,  FUN=function(i) MSEtoLogLikelihood(runif(30,1,50), runif(30,1,50))))
ggplot(modelplottingDF, aes(x =modelName, y = 1 - (nLL/randomNLL), fill = modelName, color = modelName))+
  #geom_hline(yintercept =  19.08, color = 'black', linetype = 'dotted')+
  #geom_line(aes(group=id), alpha = 0.1, color = 'black')+
  #geom_quasirandom(aes(alpha = bestModel), position = position_jitter(width = .15), size = 1,  varwidth=T) +
  geom_quasirandom( position = position_jitter(width = .15), size = 1,  varwidth=T) +
  geom_boxplot(width = .2, guides = FALSE, outlier.shape = NA, alpha = 0, color = 'black') +
  stat_summary(fun.y = mean, geom='point', shape = 5, color = 'black', fill = NA)+
  scale_alpha(range = c(0.3, 0.9))+
  #scale_y_continuous(limits = c(0,30))+
  theme_classic()+
  ylab('Log Loss')+
  xlab('Model')+
  scale_color_manual(values =c( "#009E73", "#E69F00", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), name = '')


#calculate PXP
write.table(cbind(subset(modelplottingDF, modelName == 'GP')$nLL, subset(modelplottingDF, modelName == 'dNN')$nLL, subset(modelplottingDF, modelName == 'kNN')$nLL), 'modelResults/Exp1diffevidence.csv', row.names = FALSE, sep=',', col.names=FALSE)
#Intermediate step: run pxp.m
pxpValues <- read.csv('modelResults/Exp1PXP.csv', header = F)
#build dataframe with output
pxpDF <- data.frame(model= c('GP', 'dNN', 'kNN'),  pxp= t(pxpValues))
#Reorder levels
pxpDF$model <- factor(pxpDF$model , levels= c('GP', 'dNN', 'kNN'))

#PXP
p1alt <- ggplot(pxpDF, aes(x = model, y = pxp, fill = model))+
  geom_bar(stat="identity",color='black')+
  ylab('pxp')+
  xlab('')+
  scale_y_continuous(labels = scales::percent)+
  expand_limits(y=0)+
  theme_classic()+
  scale_fill_manual(values =c( "#009E73", "#E69F00", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), name = '')+
  scale_color_manual(values =c("#009E73", "#E69F00", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), name = '')+
  theme(legend.position="none", strip.background=element_blank(), legend.key=element_rect(color=NA))
p1alt

#Main model comparison plot
p1 <- ggplot(modelplottingDF, aes(x =modelName, y = RMSE, fill = modelName, color = modelName))+
  #geom_hline(yintercept =  19.08, color = 'black', linetype = 'dotted')+
  #geom_line(aes(group=id), alpha = 0.1, color = 'black')+
  #geom_quasirandom(aes(alpha = bestModel), position = position_jitter(width = .15), size = 1,  varwidth=T) +
  geom_quasirandom( position = position_jitter(width = .15), size = 1,  varwidth=T) +
  geom_boxplot(width = .2, guides = FALSE, outlier.shape = NA, alpha = 0, color = 'black') +
  stat_summary(fun.y = mean, geom='point', shape = 5, color = 'black', fill = NA)+
  guides(fill = FALSE) +
  guides(color = FALSE) +
  guides(alpha = FALSE) +
  scale_alpha(range = c(0.3, 0.9))+
  #scale_y_continuous(limits = c(0,30))+
  theme_classic()+
  ylab('Prediction Error (RMSE)')+
  xlab('Model')+
  scale_color_manual(values =c( "#009E73", "#E69F00", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), name = '')+
  ggtitle('Model Comparison')+
  #geom_signif(comparisons=list( c("GP", "SR"), c("GP", "dNN"), c("GP", "kNN") ), annotations=c("BF = 0.2","BF=188","BF>1000"), y_position = c(14,  17,21), col="black")+
  geom_signif(comparisons=list( c("GP", "dNN"), c("GP", "kNN"), c("dNN", "kNN") ), annotations=c("BF = 3.7 x 10^4","BF = 1.5 x 10^9", "BF = 83"), y_position = c(15, 21, 18), col="black")
  
p1


#tests comparing the different models
ttestPretty(subset(modelplottingDF, modelName=="GP")$RMSE, subset(modelplottingDF, modelName=="dNN")$RMSE, paired = T, maxBF=Inf)
ttestPretty(subset(modelplottingDF, modelName=="GP")$RMSE, subset(modelplottingDF, modelName=="kNN")$RMSE, paired = T, maxBF=Inf)
ttestPretty(subset(modelplottingDF, modelName=="dNN")$RMSE, subset(modelplottingDF, modelName=="kNN")$RMSE, paired = T, maxBF=Inf)


#Compute rank-ordered confidence
confidenceDF <- data.frame()
for (pid in unique(df$id)){
    dsub <- subset(df, id == pid)
    modelUncertainty = 'gpDFuncertainty'
    dummy <- data.frame(model = 'gpDF',id=pid, rankParticipantConfidence= rank(dsub$confidenceJudgments),rankModelUncertainty = rank(dsub[,..modelUncertainty]) )
    confidenceDF <- rbind(dummy,confidenceDF)
}
confidenceDF <- subset(confidenceDF, model == 'gpDF') #include only diffusion kernel
confidenceDF$model <- factor(confidenceDF$model)
levels(confidenceDF$model)<- 'GP'

#Rank uncertainty plot
p2<- ggplot(confidenceDF, aes(x=rankParticipantConfidence, y = rankModelUncertainty, color = model, fill= model))+
  #geom_quasirandom(varwidth = T, size = 0.5, alpha = 0.2, color = '#009E73') +
  #geom_boxplot(width = .25,  outlier.shape = NA, alpha = 0.5, color = 'black') +
  stat_summary(fun.y = mean, geom = "point", color = 'black') + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", mult = 1, color = 'black', width = 0)+
  geom_smooth(aes(group=model), fill = NA, method = 'lm',formula=y~x,  se=FALSE, fullrange=TRUE, color = '#009E73')+ 
  ylab(expression(paste("GP Uncertainty Rank")))+
  xlab('Participant Confidence Rank')+
  ggtitle('GP Uncertainty and Confidence')+
  theme_classic() + theme(strip.background = element_blank(), legend.position='none')
p2


#Mixed effects model
mgpDF<-run_model(brm(confidenceJudgments ~  gpDFuncertainty  + (1+ gpDFuncertainty| id ), data = df, cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'GPconfidence') 
summary(mgpDF)
fixef(mgpDF)
#bayes_R2(mgpDF) #  0.4298 for reversed model
#mean(posterior_samples(mgpDF, pars = "b_gpDFuncertainty", exact_match = T) < 0) #Probability of model being less than 0
#tab_model(mgpDF) #very slow
#Null
mgpDFNull<-run_model(brm(confidenceJudgments ~  1  + (1+ gpDFuncertainty| id ), data = df, cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'GPconfidenceNull') 
#bayes_factor(mgpDF, mgpDFNull) #77.05063

#generate posterior predictions
xseq <- seq(min(df$gpDFuncertainty), max(df$gpDFuncertainty), length.out = 100)
newdat <-data.frame(gpDFuncertainty = xseq)
preds <- fitted(mgpDF, re_formula = NA, newdata = newdat, probs = c(0.025, 0.975))
#create new fixed effects dataframe
fixedDF <- data.frame(gpDFuncertainty =xseq,
                      confidenceJudgments = preds[,1], lower = preds[,3], upper = preds[,4] )

df$uncertaintyBreaks<- cut(df$gpDFuncertainty, breaks= seq(0, 1, by = 0.1))
levels(df$uncertaintyBreaks) <- seq(0, 1, by = 0.1)
df$uncertaintyBreaks <- as.numeric(as.character(df$uncertaintyBreaks))

p2alt<- ggplot(df, aes(x= gpDFuncertainty, y = confidenceJudgments))+
  #geom_quasirandom(varwidth = T, size = 0.5, alpha = 0.2, color = '#009E73') +
  #geom_boxplot(width = .25,  outlier.shape = NA, alpha = 0.5, color = 'black') +
  geom_line(data = fixedDF,  size = 1, color = 'black')+ #GP is
  geom_ribbon(data = fixedDF, aes(ymin=lower, ymax = upper), fill = 'black', color = NA,  alpha = 0.2 )+
  stat_summary(aes(x=uncertaintyBreaks, color =factor(uncertaintyBreaks)),fun.y = mean, geom = "point") + 
  stat_summary(aes(x=uncertaintyBreaks, color = factor(uncertaintyBreaks)), fun.data = mean_cl_boot, geom = "errorbar", mult = 1, width = 0)+
  scale_color_viridis(discrete=TRUE, direction = -1)+
  scale_fill_viridis(discrete=TRUE, direction = -1)+
  xlab(expression(paste("GP Uncertainty (", sigma,")")))+
  ylab('Participant Confidence')+
  annotate("text", x = .5, y = 9, label = "paste(italic(b)[gpUncertainty] , \" = -1.8, 95% HPD: [-2.5, -1.1]\")", parse = TRUE)+
  ggtitle('GP Uncertainty and Confidence')+
  theme_classic() + theme(strip.background = element_blank(), legend.position='none')
p2alt

#correlate confidence with GP - uncertainty
corTestPretty(df$confidenceJudgments, sqrt(df$gpDFuncertainty))

#calculate each individual's correlation
cordf <- df %>% group_by(id) %>% summarize(gpDFcor = cor(confidenceJudgments,sqrt(gpDFuncertainty), method='kendall'))
ttestPretty(na.omit(cordf$gpDFcor), mu=0) #single sample t-test to show that correlations were meaningful

#Look at parameter estimates
parameterEstimateDF <- data.frame(id = rep(participantdf$id, length(modelNames)), modelName = rep(modelNames, each=nrow(participantdf)), 
                                  parameter = rep(c("alpha","d", "k"), each=nrow(participantdf)),
                                  participantError = rep(participantdf$participantError, length(modelNames)),  
                                  estimate = c(participantdf$alpha, participantdf$dNNdelta, participantdf$kNNkValue))
parameterEstimateDF$parameter <- factor(parameterEstimateDF$parameter, levels = c("alpha", "d", "k"))
parameterEstimateDF$modelName <- factor(parameterEstimateDF$modelName, levels = modelNames)

parameterEstimateDF <- subset(parameterEstimateDF, modelName %in% c('GP-DF', 'dNN', 'kNN')) #Subset models
parameterEstimateDF$modelName <- factor(parameterEstimateDF$modelName)
levels(parameterEstimateDF$modelName) <- c('GP', 'dNN', 'kNN')

#plot of parameter estimates
p3 <- ggplot(parameterEstimateDF, aes(x=parameter, y = estimate, color = modelName, fill = modelName))+
  geom_quasirandom(position = position_jitter(width = .1),  varwidth=T, size = 0.5, alpha = 0.8) +
  geom_boxplot(width = .2,  outlier.shape = NA, alpha = 0.5, color = 'black', fill = NA) +
  stat_summary(fun.y = mean, geom='point', shape = 5, color = 'black', fill = NA)+
  guides(fill = FALSE) +
  guides(color = FALSE) +
  facet_wrap(~modelName, scales="free")+
  #scale_y_continuous(limits = c(0,30), breaks =c(0, 5,10,15, 20))+
  theme_classic()+
  ylab('Estimate')+
  xlab('')+
  ggtitle('Parameter Estimates')+
  scale_fill_manual(values =c( "#009E73", "#E69F00", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), name = '')+
  scale_color_manual(values =c("#009E73", "#E69F00", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), name = '')+
  theme(strip.background = element_blank())
p3

#Smaller than 2
ttestPretty(subset(parameterEstimateDF, parameter == 'alpha')$estimate, mu=2, maxBF=Inf)

#Comparing parameter estimates with RMSE
p4 <- ggplot(parameterEstimateDF, aes(x=estimate, y = participantError, color = modelName, fill = modelName))+
  geom_smooth(method = 'lm')+
  geom_point(alpha = .7)+
  guides(fill = FALSE) +
  guides(color = FALSE) +
  facet_wrap(~modelName, scales="free")+
  #scale_y_continuous(limits = c(0,30), breaks =c(0, 5,10,15, 20))+
  theme_classic()+
  ylab('Participant Error (RMSE)')+
  xlab('Parameter Estimate')+
  #ggtitle('Participant Error')+
  scale_fill_manual(values =c( "#009E73", "#E69F00", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), name = '')+
  scale_color_manual(values =c("#009E73", "#E69F00", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), name = '')+
  theme(strip.background = element_blank())
p4


corTestPretty(participantdf$participantError, participantdf$alpha) #Alpha is strongly correlated with error
plot(participantdf$alpha,participantdf$participantError)

corTestPretty(participantdf$participantError, participantdf$dNNdelta) #Lower d had lower error
plot(participantdf$dNNdelta, participantdf$participantError) 

corTestPretty(participantdf$participantError, participantdf$kNNkValue,  method='kendall') #Lower k had lower error
plot(participantdf$kNNkValue, participantdf$participantError) 

############################################################
#Correspondence between predictions and ground trith 
############################################################
#Omore plot looking at the correspondence between participant predictions and the ground truth
corTestPretty(df$trueTargetValue, df$predictions, maxBF = Inf)
pGroundTruthPreds <- ggplot(df, aes(x = trueTargetValue,  y=predictions )) +
  geom_point(alpha = 0.1, color = 'black')+
  theme_classic()+
  geom_smooth(method = 'lm', color = "#e41a1c") +
  xlab('True Target Value')+
  ylab('Participant Judgments')+
  #annotate("text", x = 50, y = 12, label = "paste(italic(b)[prevReward] , \" = -0.11\", \", \", italic(BF)[10] == 2.6 %*%10^43 )", parse = TRUE)+
  annotate("text", x = 25, y = 50, label = "paste(italic(r) , \" = .59\" )", parse = TRUE)
pGroundTruthPreds

corTestPretty(df$gpDFpred, df$predictions)
pGPPreds <- ggplot(df, aes(x = predictions,  y=gpDFpred )) +
  geom_point(alpha = 0.1, color = '#009E73')+
  theme_classic()+
  geom_smooth(method = 'lm', color = "#e41a1c") +
  xlab('Participant Judgments')+
  ylab('GP Predictions')+
  #annotate("text", x = 50, y = 12, label = "paste(italic(b)[prevReward] , \" = -0.11\", \", \", italic(BF)[10] == 2.6 %*%10^43 )", parse = TRUE)+
  annotate("text", x = 25, y = 50, label = "paste(italic(r) , \" = .71\")", parse = TRUE)
pGPPreds


corTestPretty(df$dNNpred, df$predictions, maxBF = Inf)
pdNNPreds <- ggplot(df, aes(x = predictions,  y= dNNpred)) +
  geom_point(alpha = 0.1, color ="#E69F00")+
  theme_classic()+
  geom_smooth(method = 'lm', color = "#e41a1c") +
  xlab('Participant Judgments')+
  ylab('dNN Predictions')+
  #annotate("text", x = 50, y = 12, label = "paste(italic(b)[prevReward] , \" = -0.11\", \", \", italic(BF)[10] == 2.6 %*%10^43 )", parse = TRUE)+
  annotate("text", x = 25, y = 50, label = "paste(italic(r) , \" = .68\")", parse = TRUE)
pdNNPreds


corTestPretty(df$kNNpred, df$predictions, maxBF = Inf)
pkNNPreds <- ggplot(df, aes(x = predictions,  y= kNNpred)) +
  geom_point(alpha = 0.1, color ="#56B4E9")+
  theme_classic()+
  geom_smooth(method = 'lm', color = "#e41a1c") +
  xlab('Participant Judgments')+
  ylab('kNN Predictions')+
  #annotate("text", x = 50, y = 12, label = "paste(italic(b)[prevReward] , \" = -0.11\", \", \", italic(BF)[10] == 2.6 %*%10^43 )", parse = TRUE)+
  annotate("text", x = 25, y = 50, label = "paste(italic(r) , \" = .67\")", parse = TRUE)
pkNNPreds


predictionScatter <- cowplot::plot_grid(pGroundTruthPreds, pGPPreds, pdNNPreds, pkNNPreds, nrow = 1, labels = 'auto')
predictionScatter
ggsave('plots/predictionScatter.pdf', predictionScatter, width = 12, height = 3, units = 'in')


############################################################
#Put it together
############################################################
pModels<- cowplot::plot_grid(p1alt, pBestFit, p3, p2alt, nrow = 1, rel_widths = c(1,1,2,2),labels = c('d', 'e', 'f', 'g'))
pModels
ggsave(filename = 'plots/exp1bottomrow.pdf', pModels, width = 12, height = 3, useDingbats=T)


complete <-cowplot::plot_grid(pBehavior,pModels,ncol = 1, labels = c('', ''))
complete
ggsave(filename = 'plots/exp1Results.pdf', complete, width = 12, height = 6, useDingbats=T)

