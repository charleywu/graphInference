#Utilities
#Charley Wu
library(data.table)

####################################################################################
# I/O
####################################################################################
#read all rds files in a folder and return them as a data frame
readRDSfolder <- function(targetFolder){
  files <- list.files(targetFolder)
  dat <- vector("list", length(files))
  for (f in 1:length(files)){
    dat[[f]] <- readRDS(paste0(targetFolder,files[f]))
  }
  return(rbindlist(dat))
}

dataImportExp1 <- function(dataFile, normalize=TRUE, delim = ','){
  #Load data
  rawdata <- read.csv(dataFile, sep=delim)
  #remove empty rows
  rawdata <- rawdata[!grepl("NULL",rawdata$task_end) & !grepl("NULL",rawdata$task_end),] 
  #extract sampleSize, records, and dat (deeply nested json data)
  sampleSize <-nrow(rawdata)
  df <- data.frame()
  #Loop through data and build df
  for (i in 1:sampleSize){
    dat = rawdata[i,]
    #general data
    id = dat$id
    age = dat$age #todo: case if the two don't match
    gender = dat$gender
    bonus = dat$reward
    #Start and end times
    start <-  strptime(as.character(dat$task_start),"%d.%m.%y %H:%M",tz="GMT")
    end <- strptime(as.character(dat$task_end),"%d.%m.%y %H:%M",tz="GMT")
    duration <- end - start
    ##Grid data
    experimentData <- fromJSON(as.character(dat$experimentData))
    #build dataframe using simple parts of the data
    dummy <-data.frame(id = rep(id,30), age = rep(age,30), gender = rep(gender, 30), bonus = rep(bonus, 30), duration = rep(duration, 30),
               graphOrder = experimentData$graphOrder[11:40], observationNum = experimentData$observationNum,  targetNodes = experimentData$targetNodes, 
               trueTargetValue = experimentData$trueTargetValue, predictions = experimentData$predictions, minList = experimentData$minList[11:40], maxList = experimentData$maxList[11:40],
               confidenceJudgments = experimentData$confidenceJudgments, errorArray = experimentData$errorArray, round = seq(1:30))
    #observedNodes has variable length, so is a bit tricky
    observedNodes <- experimentData$observedNodes
    l <- lapply(observedNodes, function(v) { #use a binary vector to indicate which nodes were observed
      default <- rep(F, 9)
      default[v]<- T
      return(default)})
    dummy[, sapply(1:9, function(i) paste0('observation', i))] <-  do.call(rbind, l) #add 7 observation columns
    df<-rbind(df, dummy)
  }
  #summary(df)
  return(df)
}

dataImportBandit <- function(dataFile, delim = ','){
  genderValues = c('Female', 'Male', 'Non-binary/Other') # lookup table for gender
  #Load data
  rawdata <- read.csv(dataFile, sep=delim, header = T)
  #remove empty rows
  rawdata <- rawdata[!grepl("NULL",rawdata$task_end) & !grepl("NULL",rawdata$task_end),] 
  #extract sampleSize, records, and dat (deeply nested json data)
  sampleSize <-1:nrow(rawdata)
  cleanedList <- sampleSize[!sampleSize %in% c(45, 92)] #Missing data
  df <- data.frame()
  #Loop through data and build df
  for (i in cleanedList){
    dat = rawdata[i,]
    #general data
    id = dat$id
    age = dat$age #todo: case if the two don't match
    gender = genderValues[as.numeric(dat$gender)]
    bonus = dat$reward
    #Start and end times
    start <-  strptime(as.character(dat$task_start),"%Y-%m-%d %H:%M:%S",tz="GMT")
    end <- strptime(as.character(dat$task_end),"%Y-%m-%d %H:%M:%S",tz="GMT")
    duration <- end - start
    ##Grid data
    experimentData <- fromJSON(as.character(dat$experimentData), simplifyVector=T, simplifyMatrix=T)
    #RTs
    RTs <- experimentData$tscollect[,2:26] - experimentData$tscollect[,1:25]
    RTs <- cbind(rep(NA,10), RTs)
    #build dataframe using simple parts of the data
    dummy <-data.frame(id = rep(id, 260), age = rep(age, 260), gender = rep(gender,260), bonus = rep(bonus,260), duration = rep(duration,260), comprehensionQuestionTries = rep(experimentData$comprehensionQuestionTries,260),
                       trial = rep(seq(0,25),10), round = rep(seq(1:10), each = 26),
                       graphOrder = rep(experimentData$envOrder[6:15], each=26)+1, #first 4 were examples, so we remove them; the 5th graph was skipped; Add +1 to shift from base_0 to base_1
                       nodeId = c(t(experimentData$xcollect)),rewardObjective =  c(t(experimentData$zcollect)),rewardSubjective =  c(t(experimentData$zscaledcollect)),
                       prevReward = c(t(cbind(rep(NA,10),experimentData$zcollect[,1:25]))),
                       RT =  c(t(RTs)), minValue = rep(experimentData$minList[5:14], each=26), maxValue = rep(experimentData$maxList[5:14], each=26))
    df<-rbind(df, dummy)
  }
  #summary(df)
  return(df)
}

####################################################################################
# Vector opterations
####################################################################################

#normalize vector
normalizeVec <- function(x) {x / ((sum(x^2))^0.5)}

normalizeMatrix <- function(m) {m / colSums(m)}

#Normalize values to min max range
normalizeMinMax <- function(x, min, max){
  normalized = (max - min) * ((x - min(x))/(max(x) - min(x))) + min
  return(normalized)
}

#rescale values to new min max from old min max range
rescaleMinMax <- function(x, newmin, newmax, oldmin = 1, oldmax = 50){
  return((newmax-newmin)/(oldmax-oldmin)*(x-oldmax)+newmax)
}

#replace all NA values in a with b
replaceAll  <- function(a,b) { 
  a[is.na(a)] <- b
}

#constrain values of vector to min max values
floorCeiling <- function(value, minValue=1, maxValue=50){
  sapply(value, FUN=function(i) min(max(i,minValue),maxValue))
}



print.plotlist<-function(xx, layout=matrix(1:length(xx), nrow=1), more=F) {
  lty<-NULL
  if ( is.matrix(layout) ) {
    lyt <- layout
    col.widths <- rep.int(1, ncol(lyt))
    row.heights <- rep.int(1, nrow(lyt))
  } else if ( is.list(layout) ) {
    stopifnot(class(layout[[1]]) == "matrix")
    lyt <- layout[[1]]
    col.widths <- if (!is.null(layout$widths)) layout$widths else rep.int(1, ncol(lyt))
    row.heights <- if (!is.null(layout$heights)) layout$heights else rep.int(1, nrow(lyt))
  }
  stopifnot(length(col.widths)==ncol(lty))
  stopifnot(length(row.heights)==nrow(lty))
  maxi <- max(lyt)
  col.pts <- cumsum(c(0, col.widths))/sum(col.widths)
  row.pts <- rev(cumsum(c(0, rev(row.heights)))/sum(row.heights))
  for(i in 1:length(xx)) {
    j <-((i-1)%%maxi)+1
    wch <- which(lyt==j, arr.ind=T)
    stopifnot(nrow(wch)>0)
    pos <- apply(wch,2,range)
    ps <- c(col.pts[pos[1,2]], row.pts[pos[2,1]+1], col.pts[pos[2,2]+1],row.pts[pos[1,1]])
    print(
      xx[[i]], 
      position = ps,
      #split=c(rev(which(lyt==j, arr.ind=T)),rev(dim(lyt))),
      more=ifelse(j != maxi & i<length(xx), T, more)
    )
  }
  invisible(F)
}

