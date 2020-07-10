#GP tutorial

#GeneralizationGradients
#Charley Wu
rm(list=ls())
packages <- c('tidyr', "ggplot2", 'dplyr',  "RColorBrewer", 'igraph', 'cowplot', 'jsonlite', 'viridis')
lapply(packages, require, character.only = TRUE)
source('models.R')

################################################################
#Shepard's exponential generalization function
################################################################
distance <- seq(0,5, length.out = 100)
df <- data.frame(distance = distance/5, generalization = exp(-distance))

p1 <- ggplot(df, aes(x = distance, y = generalization))+
  geom_line()+
  theme_classic()+
  scale_x_continuous(breaks = c(0))+
  #scale_y_continuous(breaks = c(0,1))+
  xlab('Psychological Distance')+
  ylab('Generalization')+
  ggtitle("Shepard's Law")
p1

################################################################
#RBF Generalization
################################################################

#RBF - In euclidean distance
reps <- 10000 #number of environments to generate for each length scale
kernelFun <- rbf #choose kernel function

#pre-compute all pairwise distances between points in a 8 x 8 matrix
grid <- as.matrix(expand.grid(1:8,1:8))
pairwiseDistances <- apply(grid, 1, function(i) apply(grid, 1, function(j) dist(rbind(i, j),method='manhattan')))
uniqueDistances <- unique(c(pairwiseDistances))

#create data frame
corData <- data.frame(lambda = numeric(), correlation = numeric(), distance = numeric())

#Loop over different length scales
for (lambda in c(0.5, 1, 2)){
  #generate samples from GP-rbf priors
  sampledEnvs <- mvrnorm(n = reps, mu = rep(0,64), Sigma = kernelFun(as.matrix(expand.grid(1:8,1:8)),as.matrix(expand.grid(1:8,1:8)),c(lambda,lambda,1,.0001)))
  #calculate the strength of correlations for each distance
  correlations <- c() #placeholder
  for (d in uniqueDistances){
    pairs <- which(pairwiseDistances==d, arr.ind=T) # select all pairs of points where the pairwise distance is equal to d
    valuesA <- sampledEnvs[,c(pairs[,1])] #rewards
    valuesB <- sampledEnvs[,c(pairs[,2])]
    correlations <- rbind(correlations, cor(c(valuesA),c(valuesB)))
  }
  newDat <- data.frame(lambda = rep(lambda,length(uniqueDistances)), correlation= correlations, distance = uniqueDistances)
  corData <- rbind(corData, newDat)
}

corData$Lambda <- factor(corData$lambda)
#save(corData, file = '../data/rbfCordata.Rdata')
#load("../data/rbfCordata.Rdata")


p2a <- ggplot(corData, aes(x=distance, y = correlation, color = Lambda, shape = Lambda, linetype = Lambda)) + 
  geom_line(size = 1) +
  #geom_point(size = 2) +
  ylab("Pearson Correlation")+xlab("Euclidean Distance")+
  theme_classic()+
  #xlim(c(0,10))+
  scale_x_continuous(breaks = round(seq(0, 6, by = 1),1), limits = c(0,6))+
  scale_color_brewer(palette="Set1", direction=-1, name = bquote(lambda))+
  scale_linetype_manual(values = c( 'dotted', 'twodash','solid'), name = bquote(lambda))+
  #scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
  theme(text = element_text(size=14,  family="sans"), legend.justification = c(1, 1), legend.position = c(1, 1))+
  theme(strip.background=element_blank(), legend.key=element_rect(color=NA), legend.title=element_text(size=14))
p2a

#Illustrate similarity as Euclidean distance 
#Draw a Gaussian 
m <- c(3, 3)
sigma <- matrix(c(1,0,0,1), nrow=2)
data.grid <- expand.grid(x = seq(1, 7, length.out=100), y = seq(1, 7, length.out=100))
q.samp <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = m, sigma = sigma))

p2b <- ggplot(q.samp, aes(x=x, y=y)) + 
  geom_raster(aes(fill = prob)) +
  #geom_contour(aes(z = prob), color = 'white', alpha= 0.2)+
  coord_fixed( ratio = 1)+
  theme_classic()+
  scale_fill_viridis(name = '', begin = 0.2)+
  #scale_fill_distiller(palette = "Spectral")+
  scale_x_continuous(expand = c(0.01, 0.01))+
  scale_y_continuous(expand = c(0.01, 0.01))+
  xlab('Feature 1')+
  ylab('Feature 2')+
  #annotate(geom='text', x = 3, y = 3, label=c("s"), parse=T, size = 8)+
  #annotate(geom='text', x = 6, y = 6, label=c("s*minute"), parse=T, size = 8)+
  geom_segment(x=3.4, xend=5.6, y=3.4, yend=5.6,arrow = arrow(length = unit(0.2, "cm")))+
  theme(legend.position='none',axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
p2b

#Simulate some data

simData <- data.frame(x= c(1,15,25,36,48),y= c(1, 22, 30, 40,37)) #generate some synthetic data

p2c <-ggplot(simData)+
  #geom_ribbon(aes(x=x, ymin = mu - (1.96 *sqrt(var)) , ymax =  mu + (1.96 *sqrt(var))), data = posterior, alpha = 0.2, fill = "#377eb8") + 
  #stat_smooth(method="lm",fullrange=TRUE, color = 'black', linetype = 'dotted')+
  #geom_line(data = sampleLines, aes(x = x, y = y, group = factor(sample)),  color = "#e41a1c", alpha  = 0.15)+
  #geom_line(data = posterior, aes(x=x,y=mu) ,color = "#377eb8", size = 1)+
  #scale_x_continuous(breaks = c(0,2,4,6,8))+
  coord_cartesian(xlim = c(0, 50), ylim = c(0, 50), expand = FALSE)+
  geom_point(aes(x=x, y=y))+
  theme_classic()+
  theme(legend.position="none")+
  theme( axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  labs(x = expression(s), y = expression(paste("f(s)", sep = '')))
p2c

#compute posterior
x_steps <- seq(1, 50)
new_data <- data_frame(observation = seq_along(x_steps), x = x_steps)
lengthscale <- 8
observationVar <- 100 #setting this heuristically, since the small number
noiseVar <- 0.001

krbf <- rbf(new_data$x,  new_data$x, theta = c(lengthscale,observationVar,noiseVar)) #compute kernel
posterior <- gpr(X.test = new_data$x, X = simData$x,  Y = simData$y, mu_0 =  mean(simData$y), k = krbf)  #compute posterior, where we first subtract mean of rewards from y to set mean to 0
posterior$mu <- posterior$mu 
posterior$x <- new_data$x

#calculate posterior covariance matrix
K_star <- rbf(new_data$x, new_data$x, c(lengthscale,observationVar,0)) #noiseVar is 0 herem because it is added afterwards
K_star_T <- rbf(simData$x, new_data$x, c(lengthscale,observationVar,0))
K_T_T <- rbf(simData$x, simData$x, c(lengthscale,observationVar,0)) + (diag(nrow(simData)) * noiseVar)
K_T_star <- rbf(new_data$x, simData$x, c(lengthscale,observationVar,0))
Sigma <-K_star  - t(K_star_T)%*%chol2inv(chol(K_T_T))%*%t(K_T_star)

n_draws = 50
simLines <- mvrnorm(n=n_draws, posterior$mu, Sigma)
sampleLines <- data.frame(x = rep(new_data$x, n_draws), y = c(t(simLines)), sample = rep(seq(1:n_draws), each=50))

p2d <- ggplot(simData)+
  geom_ribbon(aes(x=x, ymin = mu - (1.96 *sqrt(var)) , ymax =  mu + (1.96 *sqrt(var))), data = posterior, alpha = 0.2, fill = "#377eb8") + 
  #stat_smooth(method="lm",fullrange=TRUE, color = 'black', linetype = 'dotted')+
  geom_line(data = sampleLines, aes(x = x, y = y, group = factor(sample)),  color = "#e41a1c", alpha  = 0.15)+
  geom_line(data = posterior, aes(x=x,y=mu) ,color = "#377eb8", size = 1)+
  #scale_x_continuous(breaks = c(0,2,4,6,8))+
  coord_cartesian(xlim = c(0, 50), ylim = c(0, 50), expand = FALSE)+
  geom_point(aes(x=x, y=y))+
  theme_classic()+
  #ggtitle('Gaussian process regression')+
  theme(legend.position="none")+
  theme( axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  labs(x = expression(s), y = expression(paste("f(s)", sep = '')))
p2d


################################################################
#Diffusion Kernel Generalization
################################################################

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

#Simulate GP correlations
alphaList <- c(0.5,1,2)
gpSimDF <- data.frame()
replications <- 10000

for (i in 1:length(graphList)){ #loop through graphs
  g <- graphList[[i]] #select graph
  #Compute GP corelations
  for (alphaValue in alphaList){
    k <- diffusionKernel(g,alpha = alphaValue) #compute kernel matrix
    prior <- mvrnorm(replications, rep(0, length(V(g))), k) #sample from the prior
    #Group pairs of nodes by distance
    distanceTable <- distances(g) #distances between each node
    distanceTable[lower.tri(distanceTable)] <- NA #set lower symmetric triangle to NA
    uniqueDistances <- unique(c(distanceTable[!is.na(distanceTable)]))
    for (d in uniqueDistances){
      pairs <- which(distanceTable==d, arr.ind=T)
      valuesA <- prior[,pairs[,1]]#rewards
      valuesB <- prior[,pairs[,2]]
      dat <- data.frame(graphNum = i, alpha = alphaValue, distance = d, correlation= cor(c(valuesA),c(valuesB)), model = 'GP')
      gpSimDF<- rbind(gpSimDF,dat)
    }
  }
}

plotgpSimDF <- gpSimDF %>% group_by(distance, alpha) %>% summarize(Correlation = mean(correlation))
plotgpSimDF$alpha <- factor(plotgpSimDF$alpha)

myPalette <- colorRampPalette(rev(brewer.pal(4, "Spectral")), space="Lab")
p3a <- ggplot(plotgpSimDF, aes(x=distance, y = Correlation, color = alpha, linetype = alpha))+
  geom_line(size = 1) +
  #geom_point(size = 2) +
  ylab("Pearson Correlation")+xlab("Graph Distance")+
  theme_classic()+
  #xlim(c(0,10))+
  scale_x_continuous(breaks = round(seq(0, max(plotgpSimDF$distance), by = 1),1), limits = c(0,6))+
  #scale_color_manual(values = myPalette(4), name = expression(alpha))+
  scale_color_viridis_d( name = expression(alpha))+
  #scale_color_viridis_d(name = expression(beta))+
  #scale_color_brewer(palette="Set1", direction=-1)+
  scale_linetype_manual(values = c( 'dotted', 'twodash','solid'),  name = expression(alpha))+
  #scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
  theme(text = element_text(size=14,  family="sans"), legend.justification = c(1, 1), legend.position = c(1, 1))+
  theme(strip.background=element_blank(), legend.key=element_rect(color=NA), legend.title=element_text(size=14))
#ggtitle("Decay of Correlation")
p3a

#plot graph
library(GGally)
g <- graphList[[40]]
targetNode <- 4
V(g)$label <- V(g)
V(g)$distance <- subset(gpSimDF, alpha == 2& graphNum==40)$correlation[distances(g)[,targetNode]+1] #color based on distance from a node
p3b <- ggnet2(g,  color = "distance", palette='OrRd', label.size = 4, legend.position = "None", edge.color = 'black')
p3b



######################################################
# Put plots together
######################################################
p3a
p <- plot_grid( NULL, p3a, NULL, NULL, p2b, p2a, p2c, p2d, ncol = 4, labels = c('a','c' ,'e', 'g','b','d','f','h'))
p

ggsave('plots/Figure1.pdf', p, width = 12, height = 6, units = 'in')
