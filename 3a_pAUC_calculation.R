
library(pROC)
library(rlist)
library(stats)
library(foreach)
library(qvalue)
library(parallel)
library(doParallel)  

# Function adapted from R package pROC version 1.17.0.1
auc.roc <- function(specSensDf,
                    # Partial auc definition
                    partial.auc=FALSE, # false (consider total area) or numeric length 2: boundaries of the AUC to consider, between 0 and 1, or 0 and 100 if percent is TRUE
                    partial.auc.focus=c("specificity", "sensitivity"), # if partial.auc is not FALSE: do the boundaries
                    partial.auc.correct=FALSE,
                    allow.invalid.partial.auc.correct = FALSE,
                    percent = FALSE,
                    ... # unused required to allow roc passing arguments to plot or ci.
) {
  if (!identical(partial.auc, FALSE)) {
    partial.auc.focus <- match.arg(partial.auc.focus)
  }
  
  # Validate partial.auc
  if (! identical(partial.auc, FALSE) & !(is.numeric(partial.auc) && length(partial.auc)==2))
    stop("partial.auc must be either FALSE or a numeric vector of length 2")
  
  # Ensure partial.auc is sorted with partial.auc[1] >= partial.auc[2]
  partial.auc <- sort(partial.auc, decreasing=TRUE)
  # Get and sort the sensitivities and specificities
  
  specSensDf <- specSensDf[with(specSensDf, order(Specificity)), ]
  se <- specSensDf$Sensitivity
  sp <- specSensDf$Specificity
  
  # Full area if partial.auc is FALSE
  if (identical(partial.auc, FALSE)) {
    if (methods::is(roc, "smooth.roc") && ! is.null(roc$smoothing.args) && roc$smoothing.args$method == "binormal") {
      coefs <- coefficients(roc$model)
      auc <- unname(pnorm(coefs[1] / sqrt(1+coefs[2]^2)) * ifelse(percent, 100^2, 1))
    }
    else {
      diffs.x <- sp[-1] - sp[-length(sp)]
      means.vert <- (se[-1] + se[-length(se)])/2
      auc <- sum(means.vert * diffs.x)
    }
  }
  # Partial area
  else {
    if (partial.auc.focus == "sensitivity") {
      # if we focus on SE, just swap and invert x and y and the computations for SP will work
      x <- rev(se)
      y <- rev(sp)
    }
    else {
      x <- sp
      y <- se
    }
    
    # find the SEs and SPs in the interval
    x.inc <- x[x <= partial.auc[1] & x >= partial.auc[2]]
    y.inc <- y[x <= partial.auc[1] & x >= partial.auc[2]]
    # compute the AUC strictly in the interval
    diffs.x <- x.inc[-1] - x.inc[-length(x.inc)]
    means.vert <- (y.inc[-1] + y.inc[-length(y.inc)])/2
    auc <- sum(means.vert * diffs.x)
    # add the borders:
    if (length(x.inc) == 0) { # special case: the whole AUC is between 2 se/sp points. Need to interpolate from both
      diff.horiz <- partial.auc[1] - partial.auc[2]
      # determine indices
      idx.hi <- match(FALSE, x < partial.auc[1])
      idx.lo <- idx.hi - 1
      # proportions
      proportion.hi <- (x[idx.hi] - partial.auc[1]) / (x[idx.hi] - x[idx.lo])
      proportion.lo <- (partial.auc[2] - x[idx.lo]) / (x[idx.hi] - x[idx.lo])
      # interpolated y's
      y.hi <- y[idx.hi] + proportion.hi * (y[idx.lo] - y[idx.hi])
      y.lo <- y[idx.lo] - proportion.lo * (y[idx.lo] - y[idx.hi])
      # compute AUC
      mean.vert <- (y.hi + y.lo)/2
      auc <- mean.vert*diff.horiz
    } else { # if the upper limit is not exactly present in SPs, interpolate
      if (!(partial.auc[1] %in% x.inc)) {
        # find the limit indices
        idx.out <- match(FALSE, x < partial.auc[1])
        idx.in <- idx.out - 1
        # interpolate y
        proportion <- (partial.auc[1] - x[idx.out]) / (x[idx.in] - x[idx.out])
        y.interpolated <- y[idx.out] + proportion * (y[idx.in] - y[idx.out])
        # add to AUC
        auc <- auc + (partial.auc[1] - x[idx.in]) * (y[idx.in] + y.interpolated)/2
      }
      if (!(partial.auc[2] %in% x.inc)) { # if the lower limit is not exactly present in SPs, interpolate
        # find the limit indices in and out
        #idx.out <- length(x) - match(TRUE, rev(x) < partial.auc[2]) + 1
        idx.out <- match(TRUE, x > partial.auc[2]) - 1
        idx.in <- idx.out + 1
        # interpolate y
        proportion <- (x[idx.in] - partial.auc[2]) / (x[idx.in] - x[idx.out])
        y.interpolated <- y[idx.in] + proportion * (y[idx.out] - y[idx.in])
        # add to AUC
        auc <- auc + (x[idx.in] - partial.auc[2]) * (y[idx.in] + y.interpolated)/2
      }
    }
  }
  
  # In percent, we have 100*100 = 10,000 as maximum area, so we need to divide by a factor 100
  if (percent)
    auc <- auc/100
  
  # Correction according to McClish DC, 1989
  if (all(!identical(partial.auc, FALSE), partial.auc.correct)) { # only for pAUC
    min <- pROC:::roc.utils.min.partial.auc(partial.auc, percent)
    max <- pROC:::roc.utils.max.partial.auc(partial.auc, percent)
    # The correction is defined only when auc >= min
    if (!allow.invalid.partial.auc.correct && auc < min) {
      warning("Partial AUC correction not defined for ROC curves below the diagonal.")
      auc <- NA
    }
    else if (percent) {
      auc <- (100+((auc-min)*100/(max-min)))/2 # McClish formula adapted for %
    }
    else {
      auc <- (1+((auc-min)/(max-min)))/2 # original formula by McClish
    }
  }
  return(auc)
}

getSensAtpVal005 <- function(df2, nEcoli.pre=NA, nHuman.pre=NA){
  
  if (is.na(nEcoli.pre)){
    totalEcoli <- nrow(df2[which(grepl("ECOLI", row.names(df2))),])
    totalHuman <- nrow(df2[which(grepl("HUMAN", row.names(df2))),])
  } else {
    totalEcoli <- nEcoli.pre
    totalHuman <- nHuman.pre
  }
  
  TP <- nrow(df2[(df2$pValue<0.05) & grepl("ECOLI", row.names(df2)),])
  sens <- TP/totalEcoli
  
  return(sens)
}

getMetricsDependentOnpVal <- function(pVal, df2, y="Sens", nEcoli.pre=NA, nHuman.pre=NA){
  
  if (is.na(nEcoli.pre)){
    totalEcoli <- nrow(df2[which(grepl("ECOLI", row.names(df2))),])
    totalHuman <- nrow(df2[which(grepl("HUMAN", row.names(df2))),])
  } else {
    totalEcoli <- nEcoli.pre
    totalHuman <- nHuman.pre
  }
  
  TP <- nrow(df2[(df2$pValue<pVal) & grepl("ECOLI", row.names(df2)),])
  FP <- nrow(df2[(df2$pValue<pVal) & grepl("HUMAN", row.names(df2)),])
  oneMinusSpec <- FP/totalHuman
  sens <- TP/totalEcoli
  df2$True <- ifelse(grepl("ECOLI", row.names(df2)), "ECOLI", "HUMAN")
  df2$Prediction <- ifelse(df2$pValue<pVal, "ECOLI", "HUMAN")
  if (y == "Sens"){
    resultLst <- list(oneMinusSpec, sens)
  } else if (y == "Kappa"){
    kappa <- psych::cohen.kappa(x=df2[colnames(df2) %in% c("True", "Prediction")])[["kappa"]]
    resultLst <- list(oneMinusSpec, kappa)
  } else if (y == "Recall"){
    precision <-  TP/(TP+FP)
    resultLst <- list(sens, precision)
  }
  
  return(resultLst)
}

getpValCurvedf <- function(stats.df, nEcoli.pre=NA, nHuman.pre=NA) {
  
  if (is.na(nEcoli.pre)){
    totalEcoli <- nrow(stats.df[which(grepl("ECOLI", row.names(stats.df))),])
    totalHuman <- nrow(stats.df[which(grepl("HUMAN", row.names(stats.df))),])
  } else {
    totalEcoli <- nEcoli.pre
    totalHuman <- nHuman.pre
  }
  
  p.curve.lst <- lapply(unique(c(0, sort(stats.df$pValue)[c(TRUE, rep(FALSE, times=9))], 1)), function(pvalue, stats.df){
    #p.curve.lst <- lapply(unique(c(0, sort(stats.df$pValue), 1)), function(pvalue, stats.df){
    
    TP <- nrow(stats.df[(stats.df$pValue<=pvalue) & grepl("ECOLI", row.names(stats.df)),])
    FP <- nrow(stats.df[(stats.df$pValue<=pvalue) & grepl("HUMAN", row.names(stats.df)),])
    oneMinusSpec <- FP/totalHuman
    sens <- TP/totalEcoli
    list(oneMinusSpec, sens)
  }, stats.df=stats.df)
  
  p.curve.df <- data.frame(matrix(unlist(p.curve.lst), nrow = length(p.curve.lst), byrow = T))
  
  colnames(p.curve.df) <- c("1-Specificity", "Sensitivity")
  if (nrow(p.curve.df[(p.curve.df$`1-Specificity` == 0 & p.curve.df$Sensitivity == 0) == TRUE, ]) == 0) {
    p.curve.df <- rbind(c(0,0), p.curve.df)
    colnames(p.curve.df) <- c("1-Specificity", "Sensitivity")
  }
  
  p.curve.df <- p.curve.df[complete.cases(p.curve.df), ]
  return(p.curve.df)
}


getRegulatedProp <- function(pValueVec) {
  if (any(pValueVec > 1, na.rm = TRUE)){ # E.g. in the case of SAM there are sometimes  p values above 1
    warning(paste0(sum(pValueVec > 1), " p-values are above 1.They are set to 1."))
    pValueVec[pValueVec > 1] <- 1
  }

  # Remove missing values
  pValueVec[is.nan(pValueVec)] <- NA
  pValueVec[is.infinite(pValueVec)] <- NA
  pValueVec <- pValueVec[!is.na(pValueVec)]

  #  "missing or infinite values in inputs are not allowed"
  result <- NA
  try({
    qobj.comb <- qvalue::qvalue(p = pValueVec)
    result <- 1 - qobj.comb$pi0
    }, silent = TRUE)

  result
}


getPartialAUCs <- function(SensSpecDF) {
  partial.aucs <- c(.8, .9, .95)
  partial.auc.corrects <- c(TRUE, FALSE)
  auc.Settings <- data.frame(expand.grid(partial.auc=partial.aucs, partial.auc.correct=partial.auc.corrects))
  
  aucs.results <- apply(auc.Settings, 1, function(row) {
    # print(row["partial.auc"])
    partial.auc <- unlist(unname(row["partial.auc"]))
    if (partial.auc == 0){
      SensSpecDF2 <- rbind(SensSpecDF, c(1, max(SensSpecDF$Sensitivity),0))
      auc <- auc.roc(SensSpecDF, partial.auc=FALSE, partial.auc.focus="specificity", 
                   partial.auc.correct=as.logical(row["partial.auc.correct"]))
    } else {
      auc <- auc.roc(SensSpecDF, partial.auc=c(1, partial.auc), partial.auc.focus="specificity", 
                     partial.auc.correct=as.logical(row["partial.auc.correct"]))
    }

    auc
  })
  names(aucs.results) <- paste0("p.pauc_", auc.Settings$partial.auc, "_correct", as.logical(auc.Settings$partial.auc.correct))
  aucs.results
}

getPartialAUCsResults <- function(stats.df, nEcoli, nHuman) {
  sensAtpVal005 <- getSensAtpVal005(stats.df, nEcoli.pre=nEcoli, nHuman.pre=nHuman)
  p.roc.df <- getpValCurvedf(stats.df, nEcoli.pre=nEcoli, nHuman.pre=nHuman)
  p.roc.df$Specificity <- 1-p.roc.df$`1-Specificity`
  aucs.results <- getPartialAUCs(p.roc.df)
  aucs.results <- c(aucs.results, sensAtpVal005=sensAtpVal005)
  aucs.results
}

getRankROCdf <- function(stats.df, rankMode=c("FCpVal", "FC", "pVal")) {
  if (rankMode == "FCpVal") {
    stats.df$rank <- rank(stats.df$pValue) + rank(-stats.df$log2FC)  # use no absolute FC value here
  } else if (rankMode == "FC") {
    stats.df$rank <- rank(-stats.df$log2FC)
  } else if (rankMode == "pVal") {
    stats.df$rank <- rank(stats.df$pValue)
  }
  
  stats.df.rank <- stats.df[ , grep( "rank|Protein", names(stats.df)) ]
  stats.df.rank.nonNA <- stats.df.rank[!is.na(stats.df.rank$rank),]
  stats.df.rank.nonNA.sorted <- stats.df.rank.nonNA[order(stats.df.rank.nonNA$rank), ] 
  
  FC.roc.df <- data.frame(`No. of all proteins` = seq(nrow(stats.df.rank.nonNA.sorted)),
                          `No. of E.coli proteins` = Reduce("+", as.numeric(grepl("ECOLI", stats.df.rank.nonNA.sorted$Protein)), accumulate = TRUE))
  colnames(FC.roc.df) <- c("No. of all proteins", "No. of E.coli proteins")

  return(FC.roc.df)
}


getRankAUC <- function(stats.df, nEcoli, nHuman, rankMode = c("FCpVal", "FC", "pVal")) {
  FC.roc.df <- getRankROCdf(stats.df, rankMode = rankMode) 

  FC.roc.df <- rbind(FC.roc.df, c(nEcoli + nHuman, max(FC.roc.df$`No. of E.coli proteins`)))
    
  fc.auc <- tryCatch({
    fc.auc <- DescTools::AUC(FC.roc.df$`No. of all proteins`, FC.roc.df$`No. of E.coli proteins`, subdivisions=1000)/
      (nEcoli*nHuman)
  }, error = function(error_condition) {
    # message(error_condition)
    fc.auc <- NA
  })
  fc.auc
}


getAllRankAUCs <- function(stats.df, nEcoli, nHuman){
  rankAUC.FCpVal <- getRankAUC(stats.df, nEcoli, nHuman, rankMode = "FCpVal")
  rankAUC.FC <- getRankAUC(stats.df, nEcoli, nHuman, rankMode = "FC")
  rankAUC.pVal <- getRankAUC(stats.df, nEcoli, nHuman, rankMode = "pVal")
  c(rankAUC.FCpVal=rankAUC.FCpVal, rankAUC.FC=rankAUC.FC, rankAUC.pVal=rankAUC.pVal)
}

# Get from protein names for the three reference protein list from RDS files
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

getRegValsPropAndPAucs <- function(combIntersect=c("combined", "intersect", "DiaWorkflowProteins"), stats.df, dia) {
  row.names(stats.df) <- stats.df$Protein
  stats.protNames <- getStatsProteinNames(combIntersect = combIntersect, stats.df, dia)
  regpValsProp <- getRegulatedProp(pValueVec=stats.protNames[["stats.df.protNames"]]$pValue) 
  stats.protNames[["stats.df.protNames"]]$pValue[is.na(stats.protNames[["stats.df.protNames"]]$pValue)] <- 1 # Replace NAs with pValue of 1
  
  paucs <- getPartialAUCsResults(stats.df = stats.protNames[["stats.df.protNames"]], 
                                      nEcoli = stats.protNames[["nEcoli"]], 
                                      nHuman = stats.protNames[["nHuman"]])
  
  rankAUCs <- getAllRankAUCs(stats.df = stats.protNames[["stats.df.protNames"]], 
                             nEcoli = stats.protNames[["nEcoli"]], 
                             nHuman = stats.protNames[["nHuman"]])
  
  res <- c(paucs, 
           regpValsProp=regpValsProp, 
           rankAUCs)
  
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
      
      res <- c(getRegValsPropAndPAucs(combIntersect = "combined", stats.df, dia=NULL) ,
               getRegValsPropAndPAucs(combIntersect = "intersect", stats.df, dia=NULL),
               getRegValsPropAndPAucs(combIntersect = "DiaWorkflowProteins", stats.df, dia=dia))
      
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
  write.csv(result.df, file=paste0("pAUC_statsdfAdditionalInfos/", sub('\\.Rdata$', '', rdata),"_statsdfAdditionalInfos_single.csv"), row.names = FALSE)
  
  session <- sessionInfo()
  sink(paste0("pAUC_statsdfAdditionalInfos/sessionInfo_pAUC_statsdfAdditionalInfos_", sub('\\.Rdata$', '', rdata), "_single.txt"))
  print(session)
  sink()
}
