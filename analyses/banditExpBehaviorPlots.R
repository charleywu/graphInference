#BanditExpBehaviorPlots.R
#analyze experimentData
rm(list=ls())

packages <- c('tidyr',"jsonlite",  'tidybayes', "ggplot2", "scales",'sjPlot',  'dplyr', 'plyr', 'brms', "RColorBrewer", 'lmerTest', 'ggridges', 'cowplot', 'brms', "igraph", 'matrixcalc', 'Matrix', 'dplyr', 'corrplot', 'pander', 'BayesFactor', 'ggbeeswarm', 'viridis', 'lsr', 'ggsignif')
invisible(lapply(packages, require, character.only = TRUE))

source('exportImportGraph.R')
source('models.R')
source('utilities.R')
source('statisticalTests.R')


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

####################################################################################
# Participant data
####################################################################################
#Set of data that needs to be recovered from experiment

df <- dataImportBandit('../data/banditTask/networkBandit.csv', delim=',')
df <- subset(df, round <10) #NOTE! Removing the bonus round from behavioral analyses
df$id <- factor(df$id)


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
# Calculate graph distance, grid distance, repeat vs. unique
####################################################################################
#define node positions for plots
gridLoc <- expand.grid(1:8,1:8)
positionScale <- 80
nodePositions <- data.frame(id = 1:prod(8), x =gridLoc$Var1*positionScale, y =gridLoc$Var2*positionScale ) 

df$graphDistance <- NA
df$spatialDistance <- NA
df$nodeDegree <-NA
for (uid in unique(df$id)){
  dat <- subset(df, id == uid)
  for (r in 1:9){
    datr <- subset(dat, round == r) #subset data for round
    g <- graphList[[datr$graphOrder[[1]]]] #load graph from that round
    distanceTable <- distanceTables[[datr$graphOrder[[1]]]] #load shortest path distances between nodes
    df[df$id == uid & df$round == r,]$graphDistance <- c(NA,distanceTable[t(rbind(datr$nodeId[2:26],datr$nodeId[1:25]))]) #graph distance
    df[df$id == uid & df$round == r,]$nodeDegree <- nodeDegrees[[datr$graphOrder[[1]]]][datr$nodeId] #node degree
    vertex_attr(g, 'x', 1)
    delta_x <-V(g)$x[datr$nodeId][2:26] - V(g)$x[datr$nodeId][1:25]
    delta_y <-V(g)$y[datr$nodeId][2:26] - V(g)$y[datr$nodeId][1:25]
    df[df$id == uid & df$round == r,]$spatialDistance <- c(NA,sqrt(delta_x^2+delta_y^2))
    #plot graph
    #V(g)$x <-  nodePositions$x-1
    #V(g)$y <-  -nodePositions$y-1
    #ggnet2(g,   mode = c("x", "y"), color = "color", palette='OrRd', label = V(g), label.size = 4, legend.position = "None", edge.color = 'black')
  }
}

#Add performance info
performanceDF <- ddply(df, ~id, plyr::summarize, performance = mean(rewardObjective))
repeatDF$performance <- performanceDF$performance
corTestPretty(repeatDF$propRepeat, repeatDF$performance, maxBF = 100)

#Repeat as a function of previous reward trials
df <- df %>% group_by(id, round) %>% mutate(prevChoice = dplyr::lag(nodeId, n = 1, default = NA)) #add previous choice
df$repeatChoice <- ifelse(df$nodeId==df$prevChoice, TRUE, FALSE) #Compute repeat choice

#Mixed effects logistic regression model
repeatPrevReward <- run_model(brm(repeatChoice ~ rewardObjective + (1 + rewardObjective|id), data =subset(df, !is.na(df$repeatChoice) & trial>0),  family = "bernoulli",cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'repeatPrevReward')
repeatPrevRewardNull <- run_model(brm(repeatChoice ~ 1 + (1 + rewardObjective|id), data =subset(df, !is.na(df$repeatChoice) & trial>0),  family = "bernoulli", cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'repeatPrevRewardNull')
#bayes_factor(repeatPrevReward,repeatPrevRewardNull) #31775620410185782409221263306164916977664.00000
#bayes_R2(repeatPrevReward)  #0.5816209
fixedTerms <- fixef(repeatPrevReward)#Look at fixed terms

#Now generate predictions, removing id as a random effect
xseq <- seq(1,100)
newdat <-data.frame(rewardObjective = xseq)
preds <- fitted(repeatPrevReward, re_formula = NA, newdata = newdat, probs = c(0.025, 0.975))
#create new fixed effects dataframe
fixedDF <- data.frame( rewardObjective =xseq, repeatChoice = preds[,1], lower = preds[,3], upper = preds[,4] )
logOdds <- repeatPrevReward %>% gather_draws(b_Intercept, b_rewardObjective ) %>% median_hdi()
print(exp(logOdds[c(".value",".lower",".upper")])) #Odds ratio
#To have plottable raw data
df$rewardBreaks<- cut(df$rewardObjective, breaks=seq(0,100))
df$rewardBreaks <- as.numeric(df$rewardBreaks)

pRepeatReward <- ggplot(df, aes(rewardObjective, repeatChoice))+
  stat_summary(data = subset(df, !is.na(rewardBreaks)), aes(x = rewardBreaks, y = as.numeric(repeatChoice)), fun.y=mean, geom='point', color = "#377eb8", alpha = 0.8)+
  geom_ribbon(data = fixedDF, aes(ymin=lower, ymax = upper), color = NA, alpha = 0.4, fill = '#e41a1c')+
  geom_line(data=fixedDF, color = '#e41a1c')+
  theme_classic()+
  xlab('Reward Value')+
  ylab('P(Repeat)')+
  annotate("text", x = 50, y = 0.9, label = "OR: 1.13, 95% HPD: [1.12,1.14]")
pRepeatReward


####################################################################################
#Learning curves
####################################################################################

#Learning over trials
trialPlot<-ggplot(df, aes(trial, rewardObjective))+
  geom_hline(yintercept= 51.34095, linetype = 'dashed')+
  stat_summary(aes(group = id),fun.y = mean, geom = 'line', color = '#377eb8', size = 0.5, alpha = 0.2)+
  #geom_smooth(aes(color = id), se=F, size = 0.5)+
  stat_summary(fun.y = mean, geom = 'line', size = 1, color = '#e41a1c')+
  stat_summary(fun.data = mean_cl_boot, geom = 'ribbon', alpha = 0.2, color=NA, fill = '#e41a1c') +
  #geom_smooth(color=  'black')+
  theme_classic()+
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), name="")+
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), name="")+
  ylab('Mean Performance')+
  xlab('Trial')+
  annotate("text", x = 12.5, y = 100, label = "paste(italic(r) , \" = .93\", \", \", italic(BF)[10] == 4.5  %*%  10^7 )", parse = TRUE)+
  theme(legend.position='none', strip.background=element_blank(), legend.key=element_blank(), legend.background=element_blank())
trialPlot

#Learning over trials
df$round <- as.numeric(df$round)
roundPlot<-ggplot(df, aes(round, rewardObjective))+
  stat_summary(aes(group = id),fun.y = mean, geom = 'line', color = '#377eb8', size = 0.5, alpha = 0.2)+
  #geom_smooth(aes(color = id), se=F, size = 0.5)+
  stat_summary(fun.y = mean, geom = 'line', size = 1, color = '#e41a1c')+
  stat_summary(fun.data = mean_cl_boot, geom = 'ribbon', alpha = 0.5, color=NA, fill = '#e41a1c') +
  #geom_smooth()+
  theme_classic()+
  scale_x_continuous(breaks = pretty_breaks(n=5) )+
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), name="")+
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), name="")+
  ylab('Mean Performance')+
  xlab('Round')+
  annotate("text", x = 4, y = 100, label = "paste(italic(r) , \" = .49\", \", \", italic(BF)[10] == 1.1)", parse = TRUE)+
  theme(legend.position='none', strip.background=element_blank(), legend.key=element_blank(), legend.background=element_blank())
roundPlot

#correlation tests
dtrial <- ddply(df, ~trial, plyr::summarize, meanReward=mean(rewardObjective))
corTestPretty(dtrial$trial, dtrial$meanReward, maxBF=Inf)
dpart <- ddply(df, ~id, plyr::summarize, meanReward=mean(rewardObjective))
ttestPretty(dpart$meanReward, mu=51.34095, maxBF=Inf)

dround <- ddply(df, ~round, plyr::summarize, meanReward=mean(rewardObjective))
corTestPretty(as.numeric(as.character(dround$round)), dround$meanReward)

#####################################################################
# Graph Distance between Choices
######################################################################

#Calculate random baseline for graph distance
sampleSize <- 10000
randomDF <- data.frame()
for (i in 1:40){ #Loop through graphs and caluclate grpah distance
  g <- graphList[[i]]
  distanceTable <- distanceTables[[i]]
  samples <- sample(seq(1:64), size = sampleSize, replace = T)
  delta_x <-V(g)$x[samples][2:sampleSize] - V(g)$x[samples][1:(sampleSize-1)]
  delta_y <-V(g)$y[samples][2:sampleSize] - V(g)$y[samples][1:(sampleSize-1)]
  spatialDistance <- sqrt(delta_x^2 + delta_y^2)
  graphDistance <- distanceTable[t(rbind(samples[2:sampleSize],samples[1:(sampleSize-1)]))]
  dummy <- data.frame(nodeId = samples[2:sampleSize], graphDistance = graphDistance, spatialDistance = spatialDistance)
  randomDF <- rbind(randomDF, dummy)
}
randomDF <- rbind(randomDF, randomDF)
randomDF$category <- rep(c('Explorer', 'Exploiter'), each = nrow(randomDF)/2)
randomDF$category  <- factor(randomDF$category , levels =c('Explorer', 'Exploiter'))

#Graph Distance Plot
graphDistancePlot <- ggplot(df, aes(x=graphDistance)) + 
  geom_histogram(aes(y = ..density..*25), position = 'dodge', binwidth=1, color='black', fill = '#1b9e77')+
  stat_density(data = randomDF, aes(y = ..density..*25), geom="line",color='black', size = .8, bw = 1) + #TODO random baseline
  #geom_density(fill=NA, size = 0.7) +
  ylab("Choices Per Round") +
  xlab("Graph Distance") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5))+
  #scale_y_continuous(breaks = seq(0, 10, by = 2))+
  #ggtitle("Locality of Sampling") +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), name="")+
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), name="")+
  theme_classic()+
  #facet_grid(~category)+
  theme(legend.position='none', strip.background=element_blank(), legend.key=element_rect(color=NA))
graphDistancePlot

#Grid Distance plot
gridDistancePlot <- ggplot(df, aes(x=spatialDistance)) + 
  geom_histogram(aes(y = ..density..*250), position = 'dodge', color='black', binwidth = 10, fill="#d95f02")+
  stat_density(data = randomDF, aes(y = ..density..*250), geom="line",color='black', size = .8) +
  #stat_density(data = randomDF, aes(y = ..density..*20), geom="line",color='black', size = .8, bw = 1) + #TODO random baseline
  #geom_density(fill=NA, size = 0.7) +
  ylab("Choices Per Round") +
  xlab("Euclidean Distance") +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), name="")+
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), name="")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5))+
  #scale_y_continuous(breaks = seq(0, 10, by = 2))+
  #ggtitle("Locality of Sampling") +
  theme_classic()+
  #facet_grid(~category)+
  theme(legend.position='none', strip.background=element_blank(), legend.key=element_rect(color=NA))
gridDistancePlot


#Correlation between grid and graph distance. Can't really interpret thisn without normalizing the scales to each other
corTestPretty(df$graphDistance, df$spatialDistance, method = 'kendall')
ggplot(df, aes(x = graphDistance, y = spatialDistance))+
  geom_jitter(alpha = 0.1, color="#377eb8" )+
  geom_smooth( color = "#e41a1c")+
  theme_classic()+
  geom_abline(slope=1, linetype = 'dashed', color = 'black')+
  annotate("text", x = 20, y = 150, label = "paste(italic(r)[tau] , \" = .92\", \", \", italic(BF) ,\" > 100\" )", parse = TRUE)

#####################################################################
# Mixed effects regression of distance on previous reward value
######################################################################

DistancePrevReward <- run_model(brm(graphDistance ~ prevReward + (1 + prevReward|id), data =subset(df, !is.na(df$graphDistance)),  cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'DistancePrevReward')
DistancePrevRewardNull <- run_model(brm(graphDistance ~ 1 + (1 + prevReward|id), data =subset(df, !is.na(df$graphDistance)),  cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'DistancePrevRewardNull')
#bayes_factor(DistancePrevReward, DistancePrevRewardNull) # 25927565102350339841343034319129581444923392.00000
#bayes_R2(DistancePrevReward) #0.4238048
fixedTerms <- fixef(DistancePrevReward)#Look at fixed terms

#Now generate predictions, removing id as a random effect
xseq <- seq(1,100)
newdat <-data.frame(prevReward = xseq)
preds <- fitted(DistancePrevReward, re_formula = NA, newdata = newdat, probs = c(0.025, 0.975))
#create new fixed effects dataframe
fixedDF <- data.frame( prevReward =xseq, graphDistance = preds[,1], lower = preds[,3], upper = preds[,4] )


pDistanceReward <- ggplot(subset(df, !is.na(df$graphDistance)), aes(prevReward, graphDistance)) +
  stat_summary(aes(x = round(prevReward)-0.5), fun.y=mean,geom='point', alpha = 0.5, color = '#377eb8')+
  #geom_hline(yintercept = mean(randomDistanceDF$distance, na.rm=T ), size = 1, color = 'black', linetype='dashed')+ 
  #geom_line(data = randDF, aes(group=id), alpha = 0.2,  color = '#377eb8')+ #Random effects
  geom_line(data = fixedDF,  size = 1, color = "#e41a1c")+ #Fixed efects
  geom_ribbon(data = fixedDF, aes(ymin=lower, ymax = upper), color = NA, alpha = 0.4, fill = "#e41a1c" )+
  #geom_abline(slope = 1, linetype = 'dashed')+
  #coord_cartesian(xlim = c(0,100))+
  xlim(c(1,100))+
  theme_classic()+
  #facet_grid(~context, labeller = as_labeller(contextLabels) )+
  xlab("Previous Reward Value")+
  ylab("Distance Between Selections")+
  #annotate("text", x = 50, y = 12, label = "paste(italic(b)[prevReward] , \" = -0.11\", \", \", italic(BF)[10] == 2.6 %*%10^43 )", parse = TRUE)+
  annotate("text", x = 50, y = 12, label = "paste(italic(b)[prevReward] , \" = -0.11, 95% HPD: [-0.12, -0.10]\")", parse = TRUE)+
  theme(legend.position=c(0, 0), legend.justification=c(0,0), strip.background=element_blank(), legend.key=element_blank(), legend.background=element_blank())
pDistanceReward

corTestPretty(subset(df, !is.na(df$graphDistance))$prevReward, subset(df, !is.na(df$graphDistance))$graphDistance)

#####################################################################
#Node degree of each choice
######################################################################

#look at the overall graph statistics
graphStatistics <- data.frame()
maxRewardDegree <- c()
minRewardDegree <- c()
for (g in graphList){
  dummy= data.frame(reward=V(g)$reward, degree = degree(g))
  graphStatistics <- rbind(graphStatistics, dummy)
  maxRewardDegree <- c(maxRewardDegree, degree(g)[which.max(V(g)$reward)])
  minRewardDegree <- c(minRewardDegree, degree(g)[which.min(V(g)$reward)])
}

nodeRewardPlot <- ggplot(graphStatistics, aes(x = degree, y = reward))+ 
  geom_quasirandom(alpha=0.2, color="#377eb8" ) + 
  geom_smooth(method = 'lm', fullrange = T, color= "#e41a1c", se=F)+
  stat_summary(fun.data = mean_cl_boot, geom='errorbar', color = 'black', width = 0.2)+
  stat_summary(fun.y = mean, geom='point', shape = 23, size = 2, fill = 'white')+
  theme_classic()+
  xlab('Node degree')+
  ylab('Reward')
nodeRewardPlot

#How did people sample based on node degree
degreeDF <-df %>% group_by(trial) %>% summarize(averageNodeDegree = mean(nodeDegree))
corTestPretty(degreeDF$trial, degreeDF$averageNodeDegree) #correlated, participants chose less connectged nodes over time

nodeDegreeTrial <-ggplot(df, aes(trial, nodeDegree))+
  #geom_jitter(alpha = 0.01)+
  geom_hline(yintercept = mean(graphStatistics$degree), linetype = 'dashed')+
  geom_smooth(color  = '#e41a1c', fill = NA)+
  stat_summary(fun.y = mean, geom='point', color = '#377eb8')+
  stat_summary(fun.data = mean_cl_boot, geom='errorbar', color = '#377eb8')+
  theme_classic()+
  coord_cartesian(ylim=c(1,3))+
  ylab('Avg. Node Degree')+
  xlab('Trial')+
  scale_x_continuous(breaks = c(0,5,10,15,20,25))+
  annotate("text", x = 7, y = 3, label = "paste(italic(r) , \" = -.95\", \", \", italic(BF) ,\" > 100\" )", parse = TRUE)
nodeDegreeTrial

#TEst to see if participant preferentially choose higher or lower degree nodes
pdegreeDF <- df%>% group_by(id)  %>% summarize(averageNodeDegree = mean(nodeDegree))
ttestPretty(x=pdegreeDF$averageNodeDegree, mu=mean(graphStatistics$degree))


degreeDF <-df %>% group_by(round) %>% summarize(averageNodeDegree = mean(nodeDegree))
corTestPretty(degreeDF$round, degreeDF$averageNodeDegree)

nodeDegreeRound <- ggplot(df, aes(round, nodeDegree))+
  #geom_jitter(alpha = 0.01)+
  #geom_smooth(color  = '#e41a1c', fill = NA)+
  stat_summary(fun.y = mean, geom='line', color="#e41a1c" )+
  stat_summary(fun.y = mean, geom='point', color = '#377eb8')+
  stat_summary(fun.data = mean_cl_boot, geom='errorbar', color = '#377eb8')+
  theme_classic()+
  ylab('Node Degree')+
  scale_y_continuous(limits = c(1,3))+
  xlab('Round')
nodeDegreeRound


######################################################################
# Eigen Centrality of nodes
######################################################################
#Add in the eigen centrality of each node
eigenCentralityVec <- rep(NA, nrow(df))
for (i in 1:nrow(df)){
  subd <- df[i,]  
  eigenCentralityVec[i] <- eigenCentrality[[subd$graphOrder]][subd$nodeId]
}
df$eigenCentrality <- eigenCentralityVec

#Histogram of sampled EC vs. ground truth of all environments
meanEC <- mean(sapply(eigenCentrality, mean)) #mean eigen centrality across all graphs
allEnvs <- data.frame(eigenCentrality=as.vector(unlist(eigenCentrality)))

pEigenCentrality <- ggplot(df, aes(eigenCentrality))+
  geom_histogram(data = allEnvs, aes(y = stat(width*density),  fill = "Environment"),  alpha = 0.2)+
  geom_histogram(data = df, aes(y = stat(width*density),  fill = 'Participants'),  alpha = 0.2)+
  geom_vline(color = '#e41a1c', xintercept=meanEC, linetype = 'dashed')+
  geom_vline(color = '#377eb8', xintercept=mean(df$eigenCentrality), linetype = 'dashed')+
  theme_classic()+
  scale_fill_manual(values =c('#e41a1c','#377eb8'), name = '')+
  scale_color_manual(values =c('#e41a1c','#377eb8'), name = '')+
  scale_y_continuous(labels = percent_format()) +
  labs(x = 'Eigen Centrality', y = 'Density')+
  expand_limits(y = c(0, 0), x = c(0,0))+
  theme(legend.position=c(1,1), legend.justification = c(1,1))
pEigenCentrality  

#Do t-test
idDF <- ddply(df, ~id, plyr::summarize, EC = mean(eigenCentrality))
ttestPretty(idDF$EC, mu=meanEC, maxBF=Inf)
  
#Does EC change over trials?
ECtrial <- run_model(brm(eigenCentrality ~ trial + (1 + trial|id), data =subset(df, !is.na(df$eigenCentrality)),  cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'ECtrial')
ECtrialNull <- run_model(brm(eigenCentrality ~ 1 + (1 + trial|id), data =subset(df, !is.na(df$eigenCentrality)),  cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'ECtrialNull')
#bayes_factor(ECtrial, ECtrialNull) #158874583531.47815
#bayes_R2(ECtrial) #slow
fixedTerms <- fixef(ECtrial)#Look at fixed terms

#Now generate predictions, removing id as a random effect
xseq <- seq(0,25)
newdat <-data.frame(trial = xseq)
preds <- fitted(ECtrial, re_formula = NA, newdata = newdat, probs = c(0.025, 0.975))
#create new fixed effects dataframe
fixedDF <- data.frame( trial =xseq, eigenCentrality = preds[,1], lower = preds[,3], upper = preds[,4] )

trialwiseDF <- ddply(df, ~id+trial, plyr::summarize, eigenCentrality = mean(eigenCentrality))
pECtrial<- ggplot(trialwiseDF, aes(trial, eigenCentrality))+
  geom_hline(yintercept = meanEC, linetype = 'dashed', color='black')+
  geom_line(data = fixedDF,  size = 1, color = "#e41a1c")+ #Fixed efects
  geom_ribbon(data = fixedDF, aes(ymin=lower, ymax = upper), color = NA, alpha = 0.4, fill = "#e41a1c" )+
  stat_summary(fun.y = mean, geom='point', color = "#377eb8")+
  stat_summary(fun.data = mean_cl_boot, geom='errorbar', color = "#377eb8", width = 0)+
  #annotate("text", x = 13, y = 0.18, label = "paste(italic(b)[Trial] , \" = -.003\", \", \", italic(BF)[10] == 1.6 %*% 10^11)", parse = TRUE)+
  annotate("text", x = 13, y = 0.20, label = "paste(italic(b)[Trial] , \" = -0.003, 95% HPD: [-0.003, -0.003]\")", parse = TRUE)+
  labs(x='Trial', y = 'Eigen Centrality')+
  theme_classic()
pECtrial


#Ground truth relationship between eigenCentrality and reward
allEnvs$reward <- as.vector(sapply(graphList, function(g) V(g)$reward))
corTestPretty(allEnvs$eigenCentrality, allEnvs$reward)
pAllEnvsEC <- ggplot(allEnvs, aes(x = eigenCentrality, y = reward))+
  stat_summary(fun.y = mean, geom='point',  color = "#377eb8", alpha = 0.2)+
  geom_smooth(color =  "#e41a1c", fill = "#e41a1c")+
  theme_classic()+
  ylim(c(1,100))+
  annotate("text", x = 0.5, y = 95, label = "paste(italic(r) , \" = -.03\", \", \", italic(BF)[10] == .17 )", parse = TRUE)+
  labs(x= 'Eigen Centrality', y = 'Mean Reward')
pAllEnvsEC

#Average EC of max and min reward in each graph
maxNodes <-minNodes<- rep(NA, 40)

for (i in 1:length(eigenCentrality)){
  ECs <- eigenCentrality[[i]]
  g <- graphList[[i]]
  rewards <- V(g)$reward
  maxNodes[i] <- ECs[which.max(rewards)]
  minNodes[i] <- ECs[which.min(rewards)]
}

mean(maxNodes)
mean(minNodes)

#But Is there a relationship for the nodes that participants selected?
ECreward <- run_model(brm(rewardObjective ~ eigenCentrality + (1 + eigenCentrality|id), data =subset(df, !is.na(df$eigenCentrality) & trial>0),  cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'ECreward')
ECrewardNull <- run_model(brm(rewardObjective ~ 1 + (1 + eigenCentrality|id), data =subset(df, !is.na(df$eigenCentrality) & trial>0),  cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'ECrewardNull')
#bayes_factor(ECreward, ECrewardNull)  189620355854893547520.00000
#bayes_R2(ECreward) #0.1659
fixedTerms <- fixef(ECreward)#Look at fixed terms

#Now generate predictions, removing id as a random effect
xseq <- seq(0,1, length.out=100)
newdat <-data.frame(eigenCentrality = xseq)
preds <- fitted(ECreward, re_formula = NA, newdata = newdat, probs = c(0.025, 0.975))
#create new fixed effects dataframe
fixedDF <- data.frame( eigenCentrality =xseq, rewardObjective = preds[,1], lower = preds[,3], upper = preds[,4] )

pECRewards <- ggplot(df, aes(x = eigenCentrality, y = rewardObjective))+
  stat_summary(aes(x = round(eigenCentrality*100)/100), fun.y = mean, geom='point',  color = "#377eb8", alpha = 0.5)+
  #stat_summary(aes(x = round(eigenCentrality*100)/100), fun.data = mean_cl_boot, geom='errorbar',  color = "#377eb8", width = 0, alpha = 0.5)+
  geom_line(data = fixedDF,  size = 1, color = "#e41a1c")+ #Fixed efects
  geom_ribbon(data = fixedDF, aes(ymin=lower, ymax = upper), color = NA, alpha = 0.4, fill = "#e41a1c" )+
  theme_classic()+
  ylim(c(1,100))+
  #annotate("text", x = 0.5, y = 95, label = "paste(italic(b)[eigenCentrality] , \" = -26.5\", \", \", italic(BF)[10] == 1.9 %*%10^20 )", parse = TRUE)+
  annotate("text", x = 0.5, y = 95, label = "paste(italic(b)[eigenCentrality] , \" = -26.5, 95% HPD: [-31.2, -22.0]\")", parse = TRUE)+
  labs(x= 'Eigen Centrality', y = 'Mean Reward')
pECRewards

######################################################################
# Reaction times
######################################################################

summary(df$RT)
#Compute limits of RT distribution based on Tukey's outlier removal criteria
quartiles <- quantile(log10(df[df$RT>0 & !is.na(df$RT),'RT']/1000), probs=c(.25, .75))
H <- 1.5* IQR(log10(df[df$RT>0 & !is.na(df$RT),'RT']/1000))
ranges = c(-H,H)
limits <- quartiles +  ranges

p4c <-ggplot(df[df$RT>0 & !is.na(df$RT),], aes(x = prevReward, y = RT/1000))+
  geom_point( alpha = 0.1, color = "#377eb8" )+
  geom_smooth(color="#e41a1c",  fill ="#e41a1c", alpha = 0.2)+
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), name="")+
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), name="")+
  xlab("Reward") +
  ylab("RT (seconds) to next choice") +
  #scale_x_continuous( limits = c(0,limits[2]/1000), trans = 'log')+
  scale_y_log10(limits = c(0.01,15))+
  scale_x_continuous(breaks = seq(0, 100, by = 10))+
  annotation_logticks(sides='l')+
  theme_classic()+
  #facet_grid(category~.)+
  theme(legend.position='none', strip.background=element_blank(), legend.key=element_rect(color=NA))
p4c


df$trial5 <- cut(df$trial, breaks=c(0,5,10,15,20,25))
p5 <- ggplot(df[df$RT>0 & !is.na(df$RT),], aes(x=RT, y=trial5, fill = trial5))+
  geom_density_ridges()+
  xlab('RT in ms')+
  ylab('Trial')+
  scale_x_log10(limits = c(50,10000))+
  annotation_logticks(sides='b')+
  scale_color_manual(values = colfunc(5), name="")+
  scale_fill_manual(values = colfunc(5), name="")+
  theme_classic()+scale_y_discrete(expand = c(0.01, 0))+
  theme(legend.position='none', strip.background=element_blank(), legend.key=element_rect(color=NA))
p5

pRT <- ddply(df, c('id', 'trial'), summarize, RT = mean(RT, na.rm=T))
corTestPretty(log(pRT$RT), pRT$trial)


df$round <- factor(df$round)
p6 <- ggplot(df[df$RT>0 & !is.na(df$RT),], aes(x=RT, y=round, fill=round))+
  geom_density_ridges()+
  scale_color_manual(values = colfunc(9), name="")+
  scale_fill_manual(values = colfunc(9), name="")+
  xlab('RT in ms')+
  scale_x_log10(limits = c(50,10000))+
  annotation_logticks(sides='b')+
  #facet_grid(category~.)+
  theme_classic()+scale_y_discrete(expand = c(0.01, 0))+
  theme(legend.position='none', strip.background=element_blank(), legend.key=element_rect(color=NA))+
  ylab('Round')
p6

pRT <- ddply(df, c('id', 'round'), summarize, RT = mean(RT, na.rm=T))
corTestPretty(log(pRT$RT), as.numeric(pRT$round))


#Individual participant quartile splits
df$prevRewardValue <- NA
for (pid in unique(df$id)){
  subd <- subset(df, id==pid)
  xs <- quantile(subd$prevReward, probs=0:4/4, na.rm = T)
  df[df$id==pid,'prevRewardValue'] <- cut(subd$prevReward, breaks=xs,labels=c("Q1", "Q2", "Q3", "Q4") )
}

df$prevRewardValue <- factor(df$prevRewardValue)
levels(df$prevRewardValue) <- c("Q1", "Q2", "Q3", "Q4") 
p7 <- ggplot(subset(df, RT>0 & !is.na(prevRewardValue)), aes(x=RT, y=prevRewardValue, fill = prevRewardValue))+
  geom_density_ridges()+
  xlab('RT in ms')+
  scale_color_manual(values = colfunc(5), name="")+
  scale_fill_manual(values = colfunc(5), name="")+
  scale_x_log10(limits = c(50,10000))+
  annotation_logticks(sides='b')+
  theme_classic()+
  scale_y_discrete(expand = c(0.01, 0))+
  theme(legend.position='none', strip.background=element_blank(), legend.key=element_rect(color=NA))+
  ylab('Previous Reward')
p7

#Eigen Centrality and RT
df <- df %>% group_by(id, round) %>% mutate(prevEC = dplyr::lag(eigenCentrality, n = 1, default = NA)) #add previous eigenCentrality

df$prevEC5 <- cut(df$prevEC, 4)
pECRT <- ggplot(subset(df, RT>0), aes(x=RT, y=prevEC5, fill = prevEC5))+
  geom_density_ridges()+
  xlab('RT in ms')+
  scale_color_manual(values = colfunc(5), name="")+
  scale_fill_manual(values = colfunc(5), name="")+
  scale_x_log10(limits = c(50,10000))+
  annotation_logticks(sides='b')+
  theme_classic()+
  scale_y_discrete(expand = c(0.01, 0))+
  theme(legend.position='none', strip.background=element_blank(), legend.key=element_rect(color=NA))+
  ylab('Previous Reward')
pECRT

#mixed effects regression
ECrt <- run_model(brm(log(RT) ~ prevEC + trial + prevEC*trial + (1+prevEC|id), data =subset(df, !is.na(df$prevEC) & trial>0),  cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99), save_all_pars = T), modelName = 'ECrt')
ECrtNull <- run_model(brm(log(RT) ~ 1 + trial + prevEC*trial + (1+prevEC|id), data =subset(df, !is.na(df$prevEC) & trial>0),  cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99),  save_all_pars = T), modelName = 'ECrtNull')
#bayes_factor(ECrt, ECrtNull)  #0.89242
#bayes_R2(ECrt)
##bayes_R2(ECrtNull)
fixef(ECrt)
 
 ##########Demographics##############

demoDF <- df %>% group_by(id) %>% summarize(age = unique(age), gender = unique(gender), duration = unique(duration), bonus = unique(bonus))
mean(as.numeric(as.character(demoDF$age)))
sd(as.numeric(as.character(demoDF$age)))

table(demoDF$gender)

mean(demoDF$duration)
sd(demoDF$duration)

mean(demoDF$bonus)
sd(demoDF$bonus)

#Comprehension fails

comprehensionFails <- df %>% group_by(id) %>% summarize(compFails = unique(comprehensionQuestionTries) - 1)
ggplot(comprehensionFails, aes(compFails)) +geom_histogram() + theme_classic() +scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), name="")

#################################################################
# Putting it all together
#################################################################

#Full Paper plot
mainPlot <- cowplot::plot_grid(trialPlot,pRepeatReward, pDistanceReward, rows = 1, labels = 'auto')
mainPlot
ggsave('plots/Exp2MainResults.pdf', mainPlot, width = 12, height = 3, units = 'in')

#roundPlot 
SIplot <- cowplot::plot_grid(pEigenCentrality, pECtrial,pAllEnvsEC, pECRewards, labels = 'auto')
SIplot
ggsave('plots/Exp2SI.pdf', SIplot, width = 8, height = 5, units = 'in')


#Table
tab_model(DistancePrevReward, ECtrial)