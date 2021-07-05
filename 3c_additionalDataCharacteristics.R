library(stats)
# library(radiomics)
library(pcaMethods)
library(foreach)
library(qvalue)
library(rlist)
library(parallel)
library(doParallel)  
library(matrixStats) 

# calc_ functions are from radiomics package
calc_energy <- function(data){
  #TODO: Add dim check for 2D vs 3D
  return(sum(as.numeric(data)*as.numeric(data), na.rm=TRUE))
}

#' @describeIn first_order_features Entropy
#' @param base The base for which the logarithm is calculate
#' @param nbins The number of bins the histogram is discretized into
calc_entropy <- function(data, base=2, nbins=length(unique(c(data)))){
  # Break data into a hist
  im_range <- range(data, na.rm=TRUE)
  cuts <- table(cut(data, seq(im_range[1], im_range[2], by=diff(im_range)/nbins), include.lowest=TRUE))/length(data[!is.na(data)])
  
  #Logs cannot take 0 values, so let = 0 if no value
  entropy_vals <- vapply(cuts, function(data) ifelse(data != 0, data*logb(data, base=base), 0), FUN.VALUE = 1)
  return(-1*sum(entropy_vals))
}


calc_kurtosis <- function(data){
  n <- length(data[!is.na(data)])
  data <- data - mean(data, na.rm=TRUE)
  r <- n * sum(data^4, na.rm=TRUE) / (sum(data^2, na.rm=TRUE)^2)
  return(r * (1 - 1/n)^2 - 3)
}

calc_meanDeviation <- function(data){
  scale <- 1/prod(dim(data))
  mu <- mean(data, na.rm=TRUE)
  return(scale * sum(abs(data - mu), na.rm=TRUE))
}

calc_skewness <- function (data){
  data <- data[!is.na(data)]
  return(sum((data - mean(data))^3)/(length(data) * sd(data)^3))
}

calc_uniformity <- function(data, nbins=length(unique(c(data)))){
  # Break data into a hist
  data <- data[!is.na(data)]
  im_range <- range(data, na.rm=TRUE)
  cuts <- table(cut(data, seq(im_range[1], im_range[2], by=diff(im_range)/nbins), include.lowest=TRUE))/length(data)
  function_vals <- vapply(cuts, function(data) data^2, FUN.VALUE = 1)
  return(sum(function_vals))
}

calc_variance <- function(data) var(c(data), na.rm=TRUE)
calc_RMS <- function(data) sqrt(mean(data^2, na.rm=TRUE))

calculateMoreCharacteristics <- function(mtx, withNAs=TRUE){
  entropy <- calc_entropy(mtx)
  kurtosis <- calc_kurtosis(mtx)
  meanDeviation <- calc_meanDeviation(mtx)
  skewness <- calc_skewness(mtx)
  uniformity <- calc_uniformity(mtx)
  variance <- calc_variance(mtx)
  RMS <- calc_RMS(mtx)
  
  group.size <- ncol(mtx)/2

  var.groups.ratio <- median(matrixStats::rowVars(mtx[, 1:group.size], na.rm = TRUE)/matrixStats::rowVars(mtx[, (group.size+1):ncol(mtx)], na.rm = TRUE), na.rm = TRUE)
  
  if (withNAs){
    resultvec <- c(entropy,
                   kurtosis, 
                   meanDeviation,
                   skewness,
                   uniformity,
                   variance, 
                   RMS,
                   var.groups.ratio)
  } else {
    t.mtx <- t(mtx)
    t.mtx <- t.mtx[ , which(apply(t.mtx, 2, var) != 0)] # Remove zero variance columns 
    pca <- stats::prcomp(t.mtx, scale.=T)
    eigs <- pca$sdev^2
    prctPC1 <- eigs[1]/sum(eigs)
    elongation <- sqrt(eigs[2] / eigs[1]) # elongation
    flatness <- sqrt(eigs[length(eigs)]/eigs[1]) # flatness
    resultvec <- c(entropy,
                   kurtosis, 
                   meanDeviation,
                   skewness,
                   uniformity,
                   variance, 
                   RMS,
                   var.groups.ratio,
                   prctPC1, 
                   elongation,
                   flatness)
  }
  return(resultvec)
}

multiModel <- function(x, group.size){
  tmp <- data.frame(x2=unlist(as.vector(x)), y=as.factor(c(rep("X25", group.size), rep("X12", group.size))))
  m <- stats::lm(x2 ~ y, data=tmp)
  bp.res <- lmtest::bptest(m, studentize = TRUE) 
  bp.res[["p.value"]]
}

rds.names <- c("DIANN_DIANN_AI_GPF.rds", "DIANN_DIANN_AI.rds",
  "DIANN_MaxQuant.rds", "DIANN_MSFragger.rds",
  "DIANN_PROSIT_EDIA_GPF.rds", "OSW_DIANN_AI_GPF.rds",
  "OSW_MaxQuant.rds", "OSW_MSFragger.rds",
  "Skyline_DIANN_AI_GPF.rds", "Skyline_MaxQuant.rds",
  "Skyline_MSFragger.rds", "Skyline_PROSIT_EDIA_GPF.rds",
  "Spectronaut_DIANN_AI_GPF.rds", "Spectronaut_DirectDIA.rds",
  "Spectronaut_MaxQuant.rds", "Spectronaut_MSFragger.rds", "Spectronaut_PROSIT_EDIA_GPF.rds")

for (rds in rds.names)  {
  bootstrapLst <- readRDS(rds)
  
  ncores <- parallel::detectCores()-1
  cl <- makeCluster(ncores, type = "PSOCK")  
  doParallel::registerDoParallel(cl)  

  addInfoLst <- list()
  
  addInfoLst <- foreach(bootstrap = bootstrapLst) %dopar% {
    mtx <- as.matrix(bootstrap)
    mtx <- mtx[rowSums(is.na(mtx)) !=ncol(mtx), ]
    percNATotal2 <- mean(is.na(mtx)) * 100 
    
    characts.wNAs <- calculateMoreCharacteristics(mtx, withNAs=TRUE)
    names(characts.wNAs) <- c("entropy.wNAs",
                              "kurtosis.wNAs", 
                              "meanDeviation.wNAs",
                              "skewness.wNAs",
                              "uniformity.wNAs",
                              "variance.wNAs", 
                              "RMS.wNAs",
                              "var.groups.ratio.wNAs")
                              # "prctPC1.wNAs", 
                              #"elongation.wNAs",
                              #"flatness.wNAs")
    
    mtx <- mtx[rowSums(is.na(mtx)) == 0, ]
    nProteins.woNAs <- nrow(mtx) # number of proteins with no NAs
    characts.woNAs <- calculateMoreCharacteristics(mtx, withNAs=FALSE)
    names(characts.woNAs) <- c("entropy.woNAs",
                               "kurtosis.woNAs", 
                               "meanDeviation.woNAs",
                               "skewness.woNAs",
                               "uniformity.woNAs",
                               "variance.woNAs", 
                               "RMS.woNAs",
                               "var.groups.ratio.woNAs",
                               "prctPC1.woNAs", 
                               "elongation.woNAs",
                               "flatness.woNAs")
    
    group.size <- ncol(mtx)/2
    stats.lst <- apply(mtx, 1, multiModel, group.size=group.size)
    heterosc.woNAs <- sum(stats.lst < 0.05) / length(stats.lst)

    qobj <- qvalue::qvalue(p = stats.lst)
    heterosc.oneMinuspi0 <- 1 - qobj$pi0
    
    additional.characts <- c(characts.wNAs, characts.woNAs, heterosc.oneMinuspi0=heterosc.oneMinuspi0,
                             nProteins.woNAs=nProteins.woNAs, percNATotal2=percNATotal2)
    additional.characts
  }

  parallel::stopCluster(cl)  
  
  addInfoDf <- do.call(rbind.data.frame, addInfoLst)
  colnames(addInfoDf) <- names(addInfoLst[[1]])
  
  write.csv(addInfoDf, file=paste0(sub('\\.rds$', '', rds),"_additionalCharactInfos_allBootstrapDatasets.csv"))
}

session <- sessionInfo()
sink(paste0("sessionInfo_additionalCharactInfos.txt"))
print(session)
sink()