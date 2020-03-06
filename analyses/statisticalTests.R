#statisticalTests.R
#Charley Wu, 2019

library(BayesFactor)

####################################################################################
# Formatting for different test outputs
####################################################################################

#Pearson's R, Kendall's Tau, or Spearman's Rho: Trims leading zero and displays 2 significant digits
corformat <- function(val) { 
  out <- sub("^(-?)0.", "\\1.", sprintf("%.2f", val)) 
  return(out)
}

#p-values: 3 decimal places of precision; <.001 if the case
pformat <- function(val) { 
  if(val <.001){
    out <- '<.001'
  }else{
    out<-paste0('=',sub("^(-?)0.", "\\1.", sprintf("%.3f", val)))
  }
  return(out)
}


#Cohen's D to a single decimal
dFormat <- function(val) { 
  paste0('=',sprintf("%.1f", val))
}

#Bayes Factor: upper ceiling of 100 as >100, if greater than 1 then display as int, if less than zero show 2 decimal places
bfformat <- function(val, maxValue = 100){
  if (val>maxValue){
    out <- paste0('>',maxValue)
  }else if(val>10){
    out<-paste0('=',sprintf("%0.f", val))
  }else if(val>1 & val <10){
    out<-paste0('=',sprintf("%0.1f", val))
  }else{
    out<-paste0('=',sub("^(-?)0.", "\\1.", sprintf("%.2f", val)))
  }
  return(out)
}

#Beta coefficients:  Trims leading zero and displays 2 significant digits
betaformat <- function(val) { 
  out <- sub("^(-?)0.", "\\1.", sprintf("%.2f", val)) 
  return(out)
}
####################################################################################
# Statistical Tests
####################################################################################
ttestPretty <-function(x,y=NULL,mu=0, var.equal=T, paired=F, maxBF = 100){
  if (!is.null(y)){ #two sample
    freq <- t.test(x,y, var.equal=var.equal, paired = paired)
    d <- cohensD(x = x, y = y)
  }else{ #single sample
    freq <- t.test(x,mu = mu, var.equal=var.equal)
    d <- cohensD(x = x,mu = mu)
  }
  dof <- freq$parameter
  t <- sprintf("%.1f",freq$statistic) 
  p <- pformat(freq$p.value)
  if (!is.null(y)){ #single sample
    BF <-bfformat(extractBF(ttestBF(x,y, paired=paired))$bf, maxValue = maxBF)
  }else{#two sample
    BF <-bfformat(extractBF(ttestBF(x,mu = mu, paired=paired))$bf, maxValue = maxBF)
  }
  return(paste0('$t(',dof,')=',t,'$, $p',p, '$, $d',dFormat(d), '$, $BF',BF,'$'))
}


#Correlation test, for either pearson or kendall's tau
corTestPretty <-function(x,y, maxBF = 100, method = 'pearson'){
  freq <- cor.test(x,y, method = method)
  r <- corformat(freq$estimate) 
  p <- pformat(freq$p.value)
  if (method == 'pearson'){
    BF <-bfformat(extractBF(correlationBF(x,y))$bf, maxValue = maxBF)  
    output <- paste0('$r=',r,'$, $p',p,'$, $BF',BF,'$')
  }else if (method=='kendall'){
    BF <- bfformat(bfCorrieKernelKendallTau(tau = freq$estimate, n = length(x))$bf10, maxValue = maxBF)
    output <- paste0('$r_{\tau}=',r,'$, $p',p,'$, $BF',BF,'$')
  }
  return(output)
}



#################################################################################################
#####   Bayes factor for Kendall's tau, adapted from:                                     
#####   van Doorn, J.B., Ly, A., Marsman, M. & Wagenmakers, E.-J. (2018). Bayesian Inference for Kendall’s Rank Correlation Coefficient. The American Statistician, 72:4, 303-308, DOI: 10.1080/00031305.2016.1264998
#################################################################################################


# Prior specification Kendall's Tau
scaledBetaTau <- function(tau, alpha=1, beta=1){
  result <-   ((pi*2^(-2*alpha))/beta(alpha,alpha))  * cos((pi*tau)/2)^(2*alpha-1)
  return(result)
}

priorTau <- function(tau, kappa){
  scaledBetaTau(tau, alpha = (1/kappa), beta = (1/kappa))
}

priorTauPlus <- function(tau, kappa=1) {
  non.negative.index <- tau >=0
  less.than.one.index <- tau <=1
  value.index <- as.logical(non.negative.index*less.than.one.index)
  result <- tau*0
  result[value.index] <- 2*priorTau(tau[value.index], kappa)
  return(result)
}

priorTauMin <- function(tau, kappa=1) {
  negative.index <- tau <=0
  greater.than.min.one.index <- tau >= -1
  value.index <- as.logical(negative.index*greater.than.min.one.index)
  result <- tau*0
  result[value.index] <- 2*priorTau(tau[value.index], kappa)
  return(result)
}


# Posterior specification Kendall's Tau
postDensKendallTau <- function(delta,Tstar,n,kappa=1,var=var,test="two-sided"){ 
  if(test == "two-sided"){priorDens <- priorTau(delta,kappa)
  } else if(test == "positive"){priorDens <- priorTauPlus(delta,kappa)
  } else if(test == "negative"){priorDens <- priorTauMin(delta,kappa)}
  priorDens <- priorTau(delta,kappa)
  dens <- dnorm(Tstar,(1.5*delta*sqrt(n)),sd=sqrt(var))* priorDens
  return(dens)
}
posteriorTau <- function(delta,kentau,n,kappa=1,var=1,test="two-sided"){
  Tstar <- (kentau * ((n*(n-1))/2))/sqrt(n*(n-1)*(2*n+5)/18)
  var <- min(1,var)
  if(test == "two-sided"){lims <- c(-1,1)
  } else if(test == "positive"){lims <- c(0,1)
  } else if(test == "negative"){lims <- c(-1,0)}
  logicalCensor <- (delta >= lims[1] & delta <= lims[2])
  dens <- logicalCensor*postDensKendallTau(delta,Tstar,n,kappa,var,test=test)/
    integrate(function(delta){postDensKendallTau(delta,Tstar,n,kappa,var,test=test)},lims[1],lims[2])$value
} 

# Bayes factor computation Kendall's Tau
bfCorrieKernelKendallTau <- function(tau, n, kappa=1, var=1, ciValue=0.95){ 
  tempList <- list(vector())
  output <- list(n=n, r=tau, bf10=NA, bfPlus0=NA, bfMin0=NA)
  output$bf10 <- priorTau(0,kappa)/posteriorTau(0,tau,n,kappa=kappa,var=var,test="two-sided")
  output$bfPlus0 <- priorTauPlus(0,kappa)/posteriorTau(0,tau,n,kappa=kappa,var=var,test="positive")
  output$bfMin0 <- priorTauMin(0,kappa)/posteriorTau(0,tau,n,kappa=kappa,var=var,test="negative")
  return(output)
}

# Compute credible intervals kendalls tau
credibleIntervalKendallTau <- function(kentau,n,kappa=1,var=1, test="two-sided", ciValue = 0.95){
  nSeqs <- 1000
  lowCI <- (1-ciValue)/2
  upCI <- (1+ciValue)/2
  taus <- seq(-1,1,length.out = (nSeqs-1))
  densVals <- posteriorTau(taus, kentau, n, kappa = kappa, var = var, test = test)
  densVals <- cumsum((densVals[1:(nSeqs-1)]+densVals[2:nSeqs])*0.5*(taus[2]-taus[1]))
  lowerCI <- taus[which(densVals>=lowCI)[1]]
  upperCI <- taus[which(densVals>=upCI)[1]]
  median <- taus[which(densVals>=0.5)[1]]
  return(list(lowerCI = lowerCI, median = median, upperCI = upperCI))
}

sampleTausA <- function(myTau,myN,nSamples = 3e3, var = 1){
  nSeqs <- 1000
  tauSamples <- NULL
  taus <- seq(-1,1,length.out = nSeqs)
  densVals <- posteriorTau(taus, myTau, myN, var = var)
  ceiling <- max(densVals)
  lowerB <- taus[which(round(densVals,digits=6) != 0 )][1]
  upperB <- rev(taus[which(round(densVals,digits=6) != 0 )])[1]
  
  while(length(tauSamples) < nSamples){
    prop <- runif(1,lowerB,upperB)
    propDens <- posteriorTau(prop, myTau, myN, var = var)
    if(propDens > runif(1,0,ceiling)){tauSamples <- c(tauSamples,prop)}
  }
  return(tauSamples)
}
