#analyze experimentData
rm(list=ls())

packages <- c('plyr','tidyr',"jsonlite","cowplot", "ggplot2", 'sjPlot',"RColorBrewer", 'lmerTest', 'brms', "igraph", 'matrixcalc', 'Matrix', 'dplyr', 'corrplot', 'pander', 'BayesFactor', 'ggbeeswarm', 'viridis', 'lsr', 'ggsignif')
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
# Participant data
####################################################################################
#Set of data that needs to be recovered from experiment

df <- dataImportExp1('../data/exp1/experiment1Full.csv', delim=',')

df$RMSE <- sqrt((df$predictions - df$trueTargetValue)^2)
participantDF <- ddply(df, ~id+observationNum, plyr::summarize, RMSE = sqrt(mean((predictions - trueTargetValue)^2)), confidenceJudgments = mean(confidenceJudgments))

#Prediction Error
p1<-ggplot(data = df, aes(x =observationNum, y= RMSE, fill = factor(observationNum))) +
  #geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_hline(yintercept = 21.21, color = 'black', linetype = 'dotted')+
  geom_quasirandom(data = participantDF, aes(y = RMSE, color = factor(observationNum)), position = position_jitter(width = .15), alpha = 0.8, varwidth=T) +
  geom_boxplot(data = participantDF, width = .4, guides = FALSE, outlier.shape = NA,alpha = 0) +
  stat_summary(data =  participantDF, fun.y = mean, geom='point', shape = 5, color = 'black', fill = NA, size = 2)+
  scale_y_continuous(limits = c(0,27), breaks =c(0, 5,10,15, 20, 25))+
  scale_x_continuous(breaks =c(3,5,7))+
  scale_fill_manual(values =c(  "#0072B2", "#D55E00", "#CC79A7"), name = '')+
  scale_color_manual(values =c( "#0072B2", "#D55E00", "#CC79A7"), name = '')+
  ylab('Prediction Error (RMSE)')+
  xlab('Number of Observed Nodes')+
  ggtitle('Prediction Error')+
  theme_classic()+
  theme(legend.position ='none')
p1


#Mixed effects model
observationsRMSE<- run_model(brm(RMSE ~ observationNum + (1+ observationNum| id ), data = df, cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'observationsRMSE')
summary(observationsRMSE)
fixef(observationsRMSE)
#bayes_R2(observationsRMSE)
#Null Model
observationsRMSENull<- run_model(brm(RMSE ~ 1 + (1+ observationNum| id ), data = df, cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'observationsRMSENull')
#Bridge Sampling
#bayes_factor(observationsRMSE, observationsRMSENull) #10868703.22605

#Confidence Judgments
p2<-ggplot(data = df, aes(x =observationNum, y= confidenceJudgments, fill = factor(observationNum))) +
  #geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_quasirandom(data = participantDF, aes(y = confidenceJudgments, color = factor(observationNum)), position = position_jitter(width = .15), alpha = 0.8, varwidth=T) +
  geom_boxplot(data = participantDF,  width = .4, guides = FALSE, outlier.shape = NA, alpha = 0) +
  stat_summary(data =  participantDF, fun.y = mean, geom='point', shape = 5, color = 'black', fill = NA, size = 2)+
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_y_continuous(limits = c(1,11), breaks =c(1, 3,5,7,9,11))+
  scale_x_continuous(breaks =c(3,5,7))+
  scale_fill_manual(values =c(  "#0072B2", "#D55E00", "#CC79A7"), name = '')+
  scale_color_manual(values =c( "#0072B2", "#D55E00", "#CC79A7"), name = '')+
  ylab('Confidence (Likert Scale)')+
  xlab('Number of Observed Nodes')+
  ggtitle('Confidence Judgments')+
  theme_classic()
p2

#Mixed effects model
observationsConfidence<- run_model(brm(confidenceJudgments ~ observationNum + (1+observationNum| id ), cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), data = df, save_all_pars = T), modelName = 'observationsConfidence')
summary(observationsConfidence)
fixef(observationsConfidence)
#bayes_R2(observationsConfidence)
#Null Model
observationsConfidenceNull<- run_model(brm(confidenceJudgments ~ 1 + (1+observationNum| id ), cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), data = df, save_all_pars = T), modelName = 'observationsConfidenceNull')
#bridge Sampling
#bayes_factor(observationsConfidence, observationsConfidenceNull) #468094370.88342

#confidence vs. error
participantConfidenceDF <- ddply(df, ~id+confidenceJudgments, plyr::summarize, RMSE = sqrt(mean((predictions - trueTargetValue)^2)))
#Add rank confidence per participant
participantConfidenceDF$rankConfidence <- NA
for (pid in unique(df$id)){
  dsub <- subset(participantConfidenceDF, id == pid)
  participantConfidenceDF[participantConfidenceDF$id==pid,'rankConfidence'] = rank(dsub$confidenceJudgments)
}


#Mixed effects model
errorConfidence <- run_model(brm(RMSE ~ confidenceJudgments + (1+confidenceJudgments| id ), cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), data = df, save_all_pars = T), modelName = 'errorConfidence')
summary(errorConfidence)
fixef(errorConfidence)
#bayes_R2(errorConfidence)
#Null Model
errorConfidenceNull <- run_model(brm(RMSE ~ 1 + (1+confidenceJudgments| id ), cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), data = df, save_all_pars = T), modelName = 'errorConfidenceNull')
#bridge Sampling
#bayes_factor(errorConfidence, errorConfidenceNull) #468094370.88342

#generate posterior predictions
xseq <- seq(1,11)
newdat <-data.frame(confidenceJudgments = rep(xseq,2))
preds <- fitted(errorConfidence, re_formula = NA, newdata = newdat, probs = c(0.025, 0.975))
#create new fixed effects dataframe
fixedDF <- data.frame(confidenceJudgments = rep(xseq,2),
                      RMSE = preds[,1], lower = preds[,3], upper = preds[,4] )

p3 <- ggplot(participantConfidenceDF, aes(x = confidenceJudgments, y = RMSE, color = factor(confidenceJudgments)))+
  #geom_quasirandom(aes(color = factor(confidenceJudgments), fill = factor(confidenceJudgments)), position = position_jitter(width = .1),  varwidth=T, size = 0.8, alpha = 0.5) +
  #geom_boxplot(width = .25,  outlier.shape = NA, alpha = 0.5, color = 'black') +
  geom_line(data = fixedDF,  size = 1, color = 'black')+ #GP is
  geom_ribbon(data = fixedDF, aes(ymin=lower, ymax = upper), color = NA, alpha = 0.2 )+
  stat_summary(fun.y = mean, geom = "point") + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar",  width = 0)+
  #scale_color_brewer(palette = "Spectral", name="Confidence") +
  #scale_fill_brewer(palette = "Spectral", name="Confidence") +
  scale_color_viridis(discrete=TRUE, direction = -1)+
  scale_fill_viridis(discrete=TRUE, direction = -1)+
  #scale_y_continuous(limits = c(0,50))+
  scale_x_continuous(limits = c(0.5,11.5), breaks =c(1, 3,5,7,9,11))+
  ylab('Prediction Error (RMSE)')+
  xlab('Participant Confidence')+
  guides(fill = FALSE) +
  guides(color = FALSE) +
  annotate("text", x = 5.5, y = 22, label = "paste(italic(b)[confidence] , \" = -0.66, 95% HPD: [-0.83, -0.49]\")", parse = TRUE)+
  ggtitle('Confidence and Error')+
  theme_classic() 
p3




################################################################################################
#Put it together
################################################################################################
pBehavior<- cowplot::plot_grid(p1,p2,p3, nrow = 1, labels = 'auto')
pBehavior
ggsave(filename = 'plots/exp1TopRow.pdf', pBehavior, width = 12, height = 3, useDingbats=T)

