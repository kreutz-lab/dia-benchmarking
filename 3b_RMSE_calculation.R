library(rlist)
library(stats)
library(foreach)
library(parallel)
library(doParallel)

rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2,na.rm = TRUE))
}

nrmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2,na.rm = TRUE))/sd(actual, na.rm = TRUE)
}

getRMSE <- function(stats.df, metric) {
  stats.dfEcoli <- stats.df[grepl("ECOLI", stats.df$Protein), ]
  stats.dfHuman <- stats.df[grepl("HUMAN", stats.df$Protein), ]
  
  if (metric == "rmse"){
    Ecoli <- rmse(actual=stats.dfEcoli$log2FC, predicted=stats.dfEcoli$log2FCPredicted)
    Human <- rmse(actual=stats.dfHuman$log2FC, predicted=stats.dfHuman$log2FCPredicted)
    HumanAndEcoli <- rmse(actual=stats.df$log2FC, predicted=stats.df$log2FCPredicted)
    
  } else if (metric == "nrmse"){
    Ecoli <- nrmse(actual=stats.dfEcoli$log2FC, predicted=stats.dfEcoli$log2FCPredicted)
    Human <- nrmse(actual=stats.dfHuman$log2FC, predicted=stats.dfHuman$log2FCPredicted)
    HumanAndEcoli <- nrmse(actual=stats.df$log2FC, predicted=stats.df$log2FCPredicted)
  }
  list(Ecoli, Human, HumanAndEcoli)
}

getStatsProteinNames <- function(combIntersect, stats.df, dia=NULL) {
  if (combIntersect =="combined"){
    proteinNames <- readRDS("combinedProteinNames.rds")
  } else if (combIntersect =="intersect"){
    proteinNames <- readRDS("intersectProteinNames.rds")
  } else if (combIntersect =="DiaWorkflowProteins"){
    proteinNames <- readRDS(paste0(dia, "_ProteinNames.rds"))
  }
  nEcoli <- length(proteinNames[grepl("_ECOLI", proteinNames)])
  nHuman <- length(proteinNames[grepl("_HUMAN", proteinNames)])
  
  if (combIntersect %in% c("combined", "intersect")){
    intersectBool <- apply(stats.df, 1, function(x) {
      length(intersect(unlist(base::strsplit(x[1], ";")), proteinNames))>0
    })
    stats.df.protNames <- stats.df[intersectBool, ]
  }  else if (combIntersect =="DiaWorkflowProteins"){
    stats.df.protNames <- stats.df
  }
  
  list(stats.df.protNames=stats.df.protNames, nEcoli=nEcoli, nHuman=nHuman)
}


getCorrectedRMSEResults <- function(combIntersect=c("intersect", "DiaWorkflowProteins"), stats.df, dia) {
  
  stats.protNames <- getStatsProteinNames(combIntersect = combIntersect, stats.df, dia)

  stats.df2 <- stats.protNames[["stats.df.protNames"]]
  stats.df2[grepl("ECOLI", stats.df2$Protein), ]$log2FCPredicted <- log2(0.24038462/0.11076923) #  1.11778738, was previously wrongly assumed to be 1.058894 
  
  RMSEs <- getRMSE(stats.df2, metric="rmse")
  RMSEEcoli <- RMSEs[[1]]
  RMSEHuman <- RMSEs[[2]]
  RMSEHumanAndEcoli <- RMSEs[[3]]
  
  res <- c(RMSEEcoli=RMSEEcoli, 
           RMSEHuman=RMSEHuman, 
           RMSEHumanAndEcoli=RMSEHumanAndEcoli)
  
  names(res) <- paste0(names(res), ".", combIntersect)
  return(res)
}

getStatsForRdata <- function(rdata){  
  load(rdata)
  
  procs <- as.numeric(Sys.getenv("SLURM_NTASKS"))
  doParallel::registerDoParallel(cores=procs)
  
  stats.df.infos <- foreach::foreach(lst = result.list) %dopar% {
    
    stats.df.info <- lapply(1:length(lst), function(i){
      X1 <- i
      lst.inner <- lst[i]
      dia <- lst.inner[[1]][[2]][[1]]
      normalization <- lst.inner[[1]][[2]][[2]]
      sparsityReduction <- lst.inner[[1]][[2]][[3]]
      statTest <- lst.inner[[1]][[2]][[4]]
      groupSize <- lst.inner[[1]][[2]][[5]]
      
      stats.df <- lst.inner[[1]][[1]]
      
      res <- c(getCorrectedRMSEResults(combIntersect = "intersect", stats.df, dia=NULL),
               getCorrectedRMSEResults(combIntersect = "DiaWorkflowProteins", stats.df, dia=dia))
      
      c(X1=X1, 
        dia=dia,
        normalization=normalization, 
        sparsityReduction=sparsityReduction, 
        statTest=statTest, groupSize=groupSize,
        res
      )
    })
    stats.df.info
  }
  foreach::registerDoSEQ()
  
  stats.df.infos <- unlist(stats.df.infos, recursive = FALSE)
  
  result.df <- do.call(rbind.data.frame, stats.df.infos)
  colnames(result.df) <- names(stats.df.infos[[1]])
  write.csv(result.df, file=paste0("correctedRMSE/", sub('\\.Rdata$', '', rdata),"_correctedRMSE.csv"), row.names = FALSE)
  
  
  session <- sessionInfo()
  sink(paste0("correctedRMSE/sessionInfo_correctedRMSE_", sub('\\.Rdata$', '', rdata), "_single.txt"))
  print(session)
  sink()
  #}
  
}

