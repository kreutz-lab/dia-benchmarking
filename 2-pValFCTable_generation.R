library(MBQN)
library(mice)
library(rrcovNA)
library(limma)
library(qvalue)
library(genefilter)
library(samr)
# BiocManager::install("ROTS")
library(ROTS)
library(lasso2)
library(proDA)
library(parallel)
library(DescTools)
library(matrixcalc)
library(psych) 
library(doParallel)

#################################################################################
# NORMALIZATION
# "unnormalized", "TRQN", "QN", "median"
getDfModel <- function(modus, df) {
  if (modus == "unnormalized"){
    df.model <- df
  } else if (modus == "TRQN") {
    mtx <- as.matrix(df)
    df.trqn <- mbqn(mtx, FUN = mean)
    row.names(df.trqn) <- row.names(df)
    df.model <- as.data.frame(df.trqn)
  } else if (modus == "QN"){
    mtx <- as.matrix(df)
    df.qn <- mbqn(mtx, FUN = NULL)
    row.names(df.qn) <- row.names(df)
    df.model <- as.data.frame(df.qn)
  }else if  (modus == "median"){
    mtx <- as.matrix(df)
    df.median <- limma::normalizeMedianValues(mtx)
    df.model <- as.data.frame(df.median)
  } else {
    print("Undefined modus")
  }
  return(df.model)
}

#################################################################################
# IMPUTATION

impSeq2 <- function (x) {
  if (is.data.frame(x)) {
    x <- data.matrix(x)
  } else if (!is.matrix(x)) {
    x <- matrix(x, length(x), 1, dimnames = list(names(x), 
                                                 deparse(substitute(x))))
  }
  xcall <- match.call()
  n <- nrow(x)
  p <- ncol(x)
  isnanx = is.na(x) + 0
  risnanx = apply(isnanx, 1, sum)
  if (length(which(risnanx > 0)) == 0) 
    return(x)
  sortx <- sort.int(risnanx, index.return = TRUE)
  sorth <- sort.int(sortx$ix, index.return = TRUE)
  x = x[sortx$ix, ]
  isnanx = is.na(x) + 0
  risnanx = sortx$x
  complobs = which(risnanx == 0)
  misobs = which(risnanx != 0)
  nmisobs = length(misobs)
  ncomplobs = length(complobs)
  for (inn in 1:nmisobs) {
    # print(inn)
    if (inn == 1) {
      covx = cov(x[complobs, ])
      mx = colMeans(x[complobs, ])
    } else {
      mxo = mx
      mx = ((ncomplobs - 1) * mx + x[misobs[inn - 1], ])/ncomplobs
      covx = (ncomplobs - 2)/(ncomplobs - 1) * covx + 1/(ncomplobs - 
                                                           1) * as.matrix(x[misobs[inn - 1], ] - mx) %*% 
        t(as.matrix(x[misobs[inn - 1], ] - mx)) + as.matrix(mxo - 
                                                              mx) %*% t(as.matrix(mxo - mx))
    }
    if (p >= length(complobs) | matrixcalc::is.singular.matrix(covx)) { # added "| is.singular.matrix(covx)" due to problem with solve if covx is singular
      icovx = solve(covx + 0.01 * diag(p))
    } else {
      icovx = solve(covx)
    }
    mvar = as.logical(isnanx[misobs[inn], ])
    xo = x[misobs[inn], !mvar]
    x[misobs[inn], mvar] = mx[mvar] - solve(icovx[mvar, mvar]) %*% 
      icovx[mvar, !mvar] %*% as.matrix(xo - mx[!mvar])
    complobs = c(complobs, misobs[inn])
    ncomplobs = ncomplobs + 1
  }
  x[sorth$ix, ]
}

getSparcityReducedDf <- function(modus, df) {
  if (modus == "NoSR"){
    imp.df <- df
  } else if (modus == "SR66") {
    # Filtering 66%
    imp.df <- df[which(rowMeans(!is.na(df)) > 0.66), ]
  } else if (modus == "SR90") {
    # Filtering 90%
    imp.df <- df[which(rowMeans(!is.na(df)) > 0.9), ]
  } else {
    print("Undefined modus")
  }
  return(imp.df)
}
#################################################################################
# STATISTICAL TESTS

multiModel <- function(x, group.size, modelType){
  # print(x)
  tmp <- data.frame(x2=unlist(as.vector(x)), y=as.factor(c(rep("X25", group.size), rep("X12", group.size))))
  res <- tryCatch({
    if (modelType == "lm") {
      res <- as.vector(coef(summary(lm(x2 ~ y, data=tmp)))[,"Pr(>|t|)"][2])
    } else if (modelType == "glm") {
      res <- as.vector(coef(summary(glm(x2 ~ y, family  = Gamma(link = "identity"), data = tmp)))[,"Pr(>|t|)"][2])
    } else if (modelType == "lasso") {
      res <- as.vector(coef(summary(lasso2::l1ce(x2 ~ y, data=tmp)))[,"Pr(>|Z|)"][2])
    }
    # }, warning = function(warning_condition) {
    #  message(warning_condition)
  }, error = function(error_condition) {
    # message(error_condition)
    res <- NA
  })
  return(res)
}

statsLstToDf <- function(stats.lst){
  stats.lst[sapply(stats.lst, is.null)] <- NA
  stats.lst[sapply(stats.lst, is.nan)] <- NA
  #print(stats.lst)
  if(is.vector(stats.lst) & !is.list(stats.lst)) { 
    stats.df <- data.frame(Protein=names(stats.lst), pValue=stats.lst)
  } else {
    stats.df <- data.frame(Protein=names(stats.lst), pValue=do.call(rbind.data.frame, stats.lst)[,1])
  }
  return(stats.df)
}

getStatTestResultDf <- function(modus, df) {
  group.size <- ncol(df)/2
  groups <- c(rep("X25", group.size), rep("X12", group.size))
  seed <- NA
  if (modus == "limma") {
    # limma
    # library(limma)
    design <- model.matrix(~groups)
    fit <- lmFit(df, design)
    fit2 <- eBayes(fit)
    output <- topTable(fit2, coef = 2, number = nrow(df), sort.by = "none", adjust.method="BH")
    stats.df <- data.frame(Protein = row.names(df), pValue = output$P.Value)
  } else if (modus %in% c ("ttest","ttestBH") ){
    pVals <- rowttests(data.matrix(df),as.factor(groups), na.rm = T)$p.value
    if (modus == "ttest"){
      stats.df <- data.frame(Protein = row.names(df), pValue = pVals)
    } else if (modus == "ttestBH"){
      #  t-test mit B-H FDR correction
      pValsBH <- p.adjust(pVals, method = "BH", n = length(pVals))
      stats.df <- data.frame(Protein = row.names(df), pValue = pValsBH)
    }
  }else if  (modus == "ttestPerm"){
    # t-test mit Permutation-based FDR correction
    # http://jtleek.com/genstats/inst/doc/03_14_P-values-and-Multiple-Testing.html
    # library(qvalue)
    # library(genefilter)
    #set.seed(3333)
    r = 250
    groups.fac <- factor(groups)
    
    # Why does rowttests have slightly different results than t.test?
    tstats_obj = rowttests(data.matrix(df),groups.fac, na.rm = T)
    tstat0 = matrix(NA,nrow=dim(df)[1],ncol=r)
    tstat = tstats_obj$statistic
    
    for(i in 1:r){
      strain0 = sample(groups.fac)
      tstat0[,i] = rowttests(data.matrix(df),strain0, na.rm = T)$statistic
    }
    emp_pvals = empPvals(tstat,tstat0)
    stats.df <- data.frame(Protein = row.names(df), pValue = emp_pvals)
  } else if (modus == "SAM") {
    # SAM
    seed <- sample(1:1000000000, 1)
    set.seed(seed)
    groups <-  c(rep(1, group.size), rep(2, group.size))
    data <- list(x=as.matrix(df),y=groups, geneid=as.character(1:nrow(df)),
                 genenames=row.names(df), logged2=TRUE)
    invisible(capture.output(samr.obj<-samr(data, resp.type="Two class unpaired", nperms=250, random.seed = seed)))
    pv=data.frame(samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar))
    stats.df <- data.frame(Protein = row.names(pv), pValue = pv[,1])
  } else if (modus == "ROTS") {
    # ROTS
    groups <-  c(rep(0, group.size), rep(1, group.size))
    bool25 <- rowSums(!is.na(df[,1:group.size])) > 1
    bool12 <- rowSums(!is.na(df[,(group.size+1):ncol(df)])) > 1
    boolTooManyNAs <- which(bool25&bool12)
    
    seed <- sample(1:1000000000, 1)
    
    # results <- ROTS(data = df[boolTooManyNAs,], groups = groups , B = 100 , K = 500 , seed = 1234)
    results <- ROTS(data = df[boolTooManyNAs,], groups = groups , B = 100 , K = 500 , verbose = FALSE, seed=seed)
    results.df <- data.frame(results$logfc, results$pvalue) 
    stats.df <- data.frame(Protein = row.names(results.df), pValue = results.df$results.pvalue)
    stats.df <- rbind(stats.df, data.frame(Protein = row.names(df[-boolTooManyNAs,]), pValue = rep(NA, nrow(df[-boolTooManyNAs,]))))
    stats.df <- stats.df[order(match(stats.df$Protein,row.names(df))),]
    row.names(stats.df) <- NULL
  } else if (modus == "Lasso") {
    ## lasso
    seed <- sample(1:1000000000, 1)
    set.seed(seed)
    invisible(capture.output(stats.lst <- apply(df, 1, multiModel, group.size=group.size, modelType="lasso")))
    stats.df <- statsLstToDf(stats.lst)
  } else if (modus == "proDA") {
    # proDA --> BiocManager::install("proDA")
    # library(proDA)
    seed <- sample(1:1000000000, 1)
    set.seed(seed)
    stats.df <- tryCatch({
      fit <- proDA(data.matrix(df), design = as.factor(c(rep("X25", group.size), rep("X12", group.size))), verbose = FALSE)
      proDA.df <- test_diff(fit, X12 - X25)
      stats.df <- data.frame(Protein = proDA.df$name, pValue = proDA.df$pval)
    }, error = function(error_condition) {
      # message(error_condition)
      stats.df <- data.frame(Protein = row.names(df), pValue = NA)
    })
  } else if (modus == "LM") {
    ## LM
    seed <- sample(1:1000000000, 1)
    set.seed(seed)
    invisible(capture.output(stats.lst <- apply(df, 1, multiModel, group.size=group.size, modelType="lm")))
    stats.df <- statsLstToDf(stats.lst)
  } else if (modus == "GLMgamma") {
    # GLM-Gamma
    seed <- sample(1:1000000000, 1)
    set.seed(seed)
    stats.lst <- apply(df, 1, multiModel, group.size=group.size, modelType="glm")
    stats.df <- statsLstToDf(stats.lst)
    
  } else {
    print("Undefined modus")
  }
  return(list(stats.df, seed))
}

#################################################################################
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
    # Area Under Kappa (AUK) Curve
    kappa <- psych::cohen.kappa(x=df2[colnames(df2) %in% c("True", "Prediction")])[["kappa"]]
    resultLst <- list(oneMinusSpec, kappa)
  } else if (y == "Recall"){
    precision <-  TP/(TP+FP)
    resultLst <- list(sens, precision)
  }
  return(resultLst)
}

getpValCurvedf <- function(stats.df, y="Sens", nEcoli.pre=NA, nHuman.pre=NA) {
  row.names(stats.df) <- stats.df$Protein
  stats.df.pVal <- stats.df[ , grep( "pValue|Protein", names(stats.df)) ]
  if (is.na(nEcoli.pre)){
    p.curve.lst <- lapply(seq(0, 1, by = 0.005), getMetricsDependentOnpVal, df2=stats.df.pVal, y=y)
  } else {
    p.curve.lst <- lapply(seq(0, 1, by = 0.005), getMetricsDependentOnpVal, df2=stats.df.pVal, y=y, 
                          nEcoli.pre=nEcoli.pre, nHuman.pre=nHuman.pre)
  }
  
  p.curve.df <- data.frame(matrix(unlist(p.curve.lst), nrow=length(p.curve.lst), byrow=T))
  
  if (y == "Sens"){
    colnames(p.curve.df) <- c("1-Specificity", "Sensitivity")
    if (nrow(p.curve.df[(p.curve.df$`1-Specificity`==0 & p.curve.df$Sensitivity == 0)==TRUE, ]) == 0){
      p.curve.df <- rbind(c(0,0), p.curve.df)
      colnames(p.curve.df) <- c("1-Specificity", "Sensitivity")
    }
  } else if (y == "Kappa") {
    colnames(p.curve.df) <- c("1-Specificity", "CohensKappa")
  } else if (y == "Recall"){
    colnames(p.curve.df) <- c("Sensitivity", "Precision")
    p.curve.df <- rbind(c(0,1), p.curve.df)
    colnames(p.curve.df) <- c("Sensitivity", "Precision")
  } else if (y == "SensRank"){
  }
  p.curve.df <- p.curve.df[complete.cases(p.curve.df), ]
  return(p.curve.df)
}

getTPdependentOnFC <- function(top, df2){
  df2.sorted <- df2[order(-abs(df2$log2FC)), ]
  tops <-row.names(df2.sorted[1:top,])
  numEcoli <- sum(grepl("ECOLI", tops))
  return(list(top, numEcoli))
}

getTPdependentOnFCpValRanksum <- function(top, df2){
  df2.sorted <- df2[order(df2$FCpValRanksum), ] # use no absolute value here
  tops <-row.names(df2.sorted[1:top,])
  numEcoli <- sum(grepl("ECOLI", tops))
  return(list(top, numEcoli))
}

getFCROCdf <- function(stats.df) {
  row.names(stats.df) <- stats.df$Protein
  stats.df.log2FC <- stats.df[ , grep( "log2FC|Protein", names(stats.df)) ]
  stats.df.log2FC.nonNA <- stats.df.log2FC[!is.na(stats.df.log2FC$log2FC),]
  FC.roc.lst <- lapply(seq(nrow(stats.df.log2FC.nonNA)), getTPdependentOnFC, df2=stats.df.log2FC.nonNA)
  FC.roc.df <- data.frame(matrix(unlist(FC.roc.lst), nrow=length(FC.roc.lst), byrow=T))
  colnames(FC.roc.df) <- c("No. of all proteins", "No. of E.coli proteins")
  return(FC.roc.df)
}

getFCpValRanksumROCdf <- function(stats.df) {
  stats.df$FCpValRanksum <- rank(stats.df$pValue) + rank(-stats.df$log2FC)
  row.names(stats.df) <- stats.df$Protein
  stats.df.rank <- stats.df[ , grep( "FCpValRanksum|Protein", names(stats.df)) ]
  stats.df.rank.nonNA <- stats.df.rank[!is.na(stats.df.rank$FCpValRanksum),]
  FC.roc.lst <- lapply(seq(nrow(stats.df.rank.nonNA)), getTPdependentOnFCpValRanksum, df2=stats.df.rank.nonNA)
  FC.roc.df <- data.frame(matrix(unlist(FC.roc.lst), nrow=length(FC.roc.lst), byrow=T))
  colnames(FC.roc.df) <- c("No. of all proteins", "No. of E.coli proteins")
  return(FC.roc.df)
}

plotRocs <- function(p.curve.df, FC.roc.df){
  plt <- ggplot(data = p.curve.df, 
                aes(x=`1-Specificity`, y= `Sensitivity`)) +
    #ggtitle(comp) +
    geom_line()  + theme_minimal()
  ggsave(plt, filename = paste0("FDR.pdf"), width = 8, height = 6)
  
  plt <- ggplot(data = FC.roc.df, 
                aes(x=`No. of all proteins`, y= `No. of E.coli proteins`)) +
    #ggtitle(comp) +
    labs(x ="Top Log2FC proteins") +
    geom_line()  + theme_minimal()
  ggsave(plt, filename = paste0("FC.pdf"), width = 8, height = 6)
}

rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2,na.rm = TRUE))
}

nrmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2,na.rm = TRUE))/sd(actual, na.rm = TRUE)
}

# Kolmogorov-Smirnov Test, 
# Returns percentage of significant samples when KS test was applied on single samples compared to all samples combined
kolSmirTestSignProp <- function(mtx) {
  pvals.mtx <- apply(mtx, 2, function(x) stats::ks.test(x, mtx)$p.value)
  cnt <- sum(pvals.mtx<0.05)
  signProportion <- cnt/length(pvals.mtx)
  return(signProportion)
}

getAUCAUKRecall <- function(stats.df, nEcoli.pre=NA, nHuman.pre=NA){
  p.roc.df <- getpValCurvedf(stats.df, y="Sens", nEcoli.pre, nHuman.pre)
  
  p.auc <- tryCatch({
    p.auc <- DescTools::AUC(p.roc.df$`1-Specificity`, p.roc.df$Sensitivity, subdivisions=1000)
  }, error = function(error_condition) {
    p.auc <- NA
  })
  
  # - Area Under Kappa (AUK) Curve is computed by varying the value of the threshold and recording the Cohen’s kappa 
  # and false positive rates for all possible thresholds. for highly class-imbalanced, or skewed,data sets.
  # (x-axis FPR, y-axis Cohen's Kappa)
  p.Kapparoc.df <- getpValCurvedf(stats.df, y="Kappa", nEcoli.pre, nHuman.pre) 
  p.auk <- tryCatch({
    p.auk <- DescTools::AUC(p.Kapparoc.df$`1-Specificity`, p.Kapparoc.df$CohensKappa, subdivisions=1000)
  }, error = function(error_condition) {
    # message(error_condition)
    p.auk <- NA
  })
  
  # - Area Under Precision Recall Curve (x-axis Recall, y-axis Precision), more sensitive to the improvements 
  # for the positive class compared to ROC AUC.  One common scenario is a highly imbalanced dataset where the fraction of positive class, which we want to find, is small
  p.PrecRecallroc.df <- getpValCurvedf(stats.df, y="Recall", nEcoli.pre, nHuman.pre)
  p.PrecRecallauc <- tryCatch({
    p.PrecRecallauc <- DescTools::AUC(p.PrecRecallroc.df$Sensitivity, p.PrecRecallroc.df$Precision, subdivisions=1000)
  }, error = function(error_condition) {
    # message(error_condition)
    p.PrecRecallauc <- NA
  })
  
  list(p.auc, p.auk, p.PrecRecallauc)
}

getFCauc <- function(stats.df, FCpValRanksum=FALSE) {
  if(FCpValRanksum){
    FC.roc.df <- getFCpValRanksumROCdf(stats.df)
  } else{
    FC.roc.df <- getFCROCdf(stats.df) 
  }

  fc.auc <- tryCatch({
    fc.auc <- DescTools::AUC(FC.roc.df$`No. of all proteins`, FC.roc.df$`No. of E.coli proteins`, subdivisions=1000)/(max(FC.roc.df$`No. of all proteins`)* max(FC.roc.df$`No. of E.coli proteins`))
  }, error = function(error_condition) {
    # message(error_condition)
    fc.auc <- NA
  })
  fc.auc
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

runAnalysisForEachBootstrapSample <- function(df, dia, normalization, sparcityReduction, statTest,
                                              combinedProteinNames, intersectProteinNames) {
  # Remove empty rows
  df <- df[rowSums(is.na(df)) != ncol(df), ]     
  
  nEcoli.pre <- sum(grepl("ECOLI", row.names(df)))
  nHuman.pre <- sum(grepl("HUMAN", row.names(df)))  
  group.size <- ncol(df)/2
  
  medianSampleVariance <- median(apply(df, 2, var,  na.rm=TRUE))
  medianProteinVariance <- median(unname(apply(df, 1, var,  na.rm=TRUE)), na.rm = TRUE)
  KS.SignProp <- kolSmirTestSignProp(as.matrix(df))
  
  percNATotal <- mean(is.na(df)) * 100 # With rows with only NAs removed before?
  percOfRowsWithNAs <- sum(apply(df, 1, anyNA))/nrow(df) * 100
  
  sparcRed.runtime <- system.time({ 
    # Sparcity Reduction
    df.sr <- getSparcityReducedDf(sparcityReduction, df) 
  })   
  
  normalization.runtime <- system.time({ 
    # Normalization
    #print("Normalization")
    df <- getDfModel(normalization, df.sr)
  })
  
  log2FC.df <- data.frame(Protein=row.names(df), log2FC=rowMeans(df[, (group.size+1):ncol(df)], na.rm = TRUE)-rowMeans(df[, 1:group.size], na.rm = TRUE))
  
  statTest.runtime <- system.time({ 
    # print("Stats")
    stats.df.lst <- getStatTestResultDf(statTest, df) 
    stats.df <- stats.df.lst[[1]]
    seed.stat <- stats.df.lst[[2]]
  })
  
  sparcRed.runtime <- sparcRed.runtime["user.self"]
  normalization.runtime <- normalization.runtime["user.self"]
  statTest.runtime <- statTest.runtime["user.self"]
  
  # Add Bonferroni Hochberg corrected p values
  stats.df$pValueBH <- p.adjust(stats.df$pValue, method="BH")
  
  stats.df <- merge(stats.df, log2FC.df, by=c("Protein"), all.x=FALSE)
  stats.df$log2FCPredicted <- NA
  stats.df[grepl("ECOLI", stats.df$Protein), ]$log2FCPredicted <- log2((1/12)/(1/25)) # 1.058894
  stats.df[grepl("HUMAN", stats.df$Protein), ]$log2FCPredicted <- log2(1)
  
  nEcoli.comb <- length(combinedProteinNames[grepl("_ECOLI", combinedProteinNames)])
  # # [1] 2127
  nHuman.comb <- length(combinedProteinNames[grepl("_HUMAN", combinedProteinNames)])
  # # [1] 11516
  # 
  nEcoli.intersect <- length(intersectProteinNames[grepl("_ECOLI", intersectProteinNames)])
  # # [1] 745
  nHuman.intersect <- length(intersectProteinNames[grepl("_HUMAN", intersectProteinNames)])
  # # [1] 4499
  
  intersectBoolComb <- apply(stats.df, 1, function(x) {
    length(intersect(unlist(base::strsplit(x[1], ";")), combinedProteinNames))>0
  })
  stats.df.combProtNames <- stats.df[intersectBoolComb, ]

  intersectBoolIntersect <- apply(stats.df, 1, function(x) {
    length(intersect(unlist(base::strsplit(x[1], ";")), intersectProteinNames))>0
  })
  stats.df.intersectProtNames <- stats.df[intersectBoolIntersect, ]
  
  AUCAUKRecallLst <- getAUCAUKRecall(stats.df, nEcoli.pre=NA, nHuman.pre=NA)
  p.auc <- AUCAUKRecallLst[[1]]
  p.auk <- AUCAUKRecallLst[[2]]
  p.PrecRecallauc <- AUCAUKRecallLst[[3]]  
  
  AUCAUKRecallLst.pre <- getAUCAUKRecall(stats.df, nEcoli.pre=nEcoli.pre, nHuman.pre=nHuman.pre)
  p.auc.pre <- AUCAUKRecallLst.pre[[1]]
  p.auk.pre <- AUCAUKRecallLst.pre[[2]]
  p.PrecRecallauc.pre <- AUCAUKRecallLst.pre[[3]]
  
  AUCAUKRecallLst.comb <- getAUCAUKRecall(stats.df.combProtNames, nEcoli.pre=nEcoli.comb, nHuman.pre=nHuman.comb)
  p.auc.comb <- AUCAUKRecallLst.comb[[1]]
  p.auk.comb <- AUCAUKRecallLst.comb[[2]]
  p.PrecRecallauc.comb <- AUCAUKRecallLst.comb[[3]]
  
  AUCAUKRecallLst.intersect <- getAUCAUKRecall(stats.df.intersectProtNames, nEcoli.pre=nEcoli.intersect, nHuman.pre=nHuman.intersect)
  p.auc.intersect <- AUCAUKRecallLst.intersect[[1]]
  p.auk.intersect <- AUCAUKRecallLst.intersect[[2]]
  p.PrecRecallauc.intersect <- AUCAUKRecallLst.intersect[[3]]
  
  fc.auc <- getFCauc(stats.df) 
  fc.auc.intersectProtNames <- getFCauc(stats.df.intersectProtNames)
  
  # Calculate RMSE
  
  RMSEs <- getRMSE(stats.df, metric="rmse")
  RMSEEcoli <- RMSEs[[1]]
  RMSEHuman <- RMSEs[[2]]
  RMSEHumanAndEcoli <- RMSEs[[3]]
  
  NRMSEs <- getRMSE(stats.df, metric="nrmse")
  NRMSEEcoli <- NRMSEs[[1]]
  NRMSEHuman <- NRMSEs[[2]]
  NRMSEHumanAndEcoli <- NRMSEs[[3]]
  
  RMSEs.Intersect <- getRMSE(stats.df.intersectProtNames, metric="rmse")
  RMSEEcoli.Intersect <- RMSEs.Intersect[[1]]
  RMSEHuman.Intersect <- RMSEs.Intersect[[2]]
  RMSEHumanAndEcoli.Intersect <- RMSEs.Intersect[[3]]
  
  NRMSEs.Intersect <- getRMSE(stats.df.intersectProtNames, metric="nrmse")
  NRMSEEcoli.Intersect <- NRMSEs.Intersect[[1]]
  NRMSEHuman.Intersect <- NRMSEs.Intersect[[2]]
  NRMSEHumanAndEcoli.Intersect <- NRMSEs.Intersect[[3]]
  
  # AUC based on rank(stats.df$pValue) + rank(-stats.df$log2FC), x-axis "No. of all proteins", y-axis "No. of E.coli proteins"
  fc.auc.FCpValRanksum <- getFCauc(stats.df, FCpValRanksum=TRUE) 
  fc.auc.FCpValRanksum.intersectProtNames <- getFCauc(stats.df.intersectProtNames, FCpValRanksum=TRUE) 
  
  nEcoli <- nrow(stats.df[which(grepl("ECOLI", stats.df$Protein)),])
  nHuman <- nrow(stats.df[which(grepl("HUMAN", stats.df$Protein)),])
  
  summary <- list(dia, normalization, sparcityReduction, statTest, groupSize=group.size, nAllProteins=nrow(stats.df), 
                  nEcoliProteins=nEcoli, nHumanProteins=nHuman, nEcoliProteins.pre=nEcoli.pre, nHumanProteins.pre=nHuman.pre,
                  medianSampleVariance=medianSampleVariance, medianProteinVariance=medianProteinVariance, KS.SignProp=KS.SignProp, percNATotal=percNATotal, percOfRowsWithNAs=percOfRowsWithNAs,
                  sparcRed.runtime = sparcRed.runtime, normalization.runtime = normalization.runtime, statTest.runtime = statTest.runtime,
                  seed.stat=seed.stat,
                  p.auc=p.auc, p.auk=p.auk, p.PrecRecallauc=p.PrecRecallauc,
                  p.auc.pre=p.auc.pre, p.auk.pre=p.auk.pre, p.PrecRecallauc.pre=p.PrecRecallauc.pre, 
                  p.auc.comb=p.auc.comb, p.auk.comb=p.auk.comb, p.PrecRecallauc.comb=p.PrecRecallauc.comb, 
                  p.auc.intersect=p.auc.intersect, p.auk.intersect=p.auk.intersect, p.PrecRecallauc.intersect=p.PrecRecallauc.intersect, 
                  fc.auc=fc.auc, fc.auc.intersectProtNames=fc.auc.intersectProtNames,
                  fc.auc.FCpValRanksum=fc.auc.FCpValRanksum, fc.auc.FCpValRanksum.intersectProtNames=fc.auc.FCpValRanksum.intersectProtNames,
                  RMSEEcoli=RMSEEcoli, RMSEHuman=RMSEHuman, RMSEHumanAndEcoli=RMSEHumanAndEcoli,
                  NRMSEEcoli=NRMSEEcoli, NRMSEHuman=NRMSEHuman, NRMSEHumanAndEcoli=NRMSEHumanAndEcoli,
                  RMSEEcoli.Intersect=RMSEEcoli.Intersect, RMSEHuman.Intersect=RMSEHuman.Intersect, RMSEHumanAndEcoli.Intersect=RMSEHumanAndEcoli.Intersect,
                  NRMSEEcoli.Intersect=NRMSEEcoli.Intersect, NRMSEHuman.Intersect=NRMSEHuman.Intersect, NRMSEHumanAndEcoli.Intersect=NRMSEHumanAndEcoli.Intersect)
  return(list(stats.df, summary))
}

runAnalysis <- function(vector, repListNorm, repListUnnorm) {
  # print(vector)
  normalization <- as.character(vector[1])
  sparcityReduction <- as.character(vector[2])
  statTest <- as.character(vector[3])
  print(normalization)
  print(sparcityReduction)
  print(statTest)
  print("--------------------------")
  
  subresult.list <- lapply(repListNorm, runAnalysisForEachBootstrapSample,
                           normalization=normalization, sparcityReduction=sparcityReduction, statTest=statTest)
  
  return(subresult.list)
}

################################################################################

# batch <- 1 # max 32 bei batch.size 20, max 64 bei batch.size 10, max 128 bei batch.size 5
run <- function(batch, batch.size) {
  
  session <- sessionInfo()
  sink(paste0("sessionInfo_", batch, "_",batch.size, ".txt"))
  print(session)
  sink()
  
  # GET ALL PARAMETER COMBINATIONS
  dias <- c("OSW_MaxQuant", "Skyline_MaxQuant", "Skyline_PROSIT_EDIA_GPF", 
            "Spectronaut_PROSIT_EDIA_GPF", "Spectronaut_DirectDIA", "Spectronaut_DIANN_AI_GPF", 
            "Skyline_DIANN_AI_GPF", "OSW_DIANN_AI_GPF", "Spectronaut_MSFragger", 
            "DIANN_PROSIT_EDIA_GPF", "DIANN_DIANN_AI", "DIANN_MSFragger", 
            "Skyline_MSFragger", "Spectronaut_MaxQuant", "OSW_MSFragger", 
            "DIANN_MaxQuant", "DIANN_DIANN_AI_GPF")
  normalizations <- c("unnormalized", "TRQN", "QN", "median")
  sparcityReductions <- c("NoSR", "SR66", "SR90")
  statTests <- c("ttest", "limma", "LM", "GLMgamma", "Lasso", "SAM", "ROTS") # lm() is conducted in MSstats
  
  combs <- expand.grid(statTests, sparcityReductions, normalizations, dias)
  colnames(combs) <- c("statTest", "sparcityReduction", "normalization", "dia")
 
  if (batch*batch.size > nrow(combs)){
    combs <- combs[(((batch-1)*batch.size)+1):nrow(combs),]
  } else{
    combs <- combs[(((batch-1)*batch.size)+1):(batch*batch.size),]
  }
  
  combs.lst <- as.list(as.data.frame(t(combs)))
  
  procs <- as.numeric(Sys.getenv("SLURM_NTASKS"))
  registerDoParallel(cores=procs)
  
  combinedProteinNames <- readRDS("combinedProteinNames.rds")
  intersectProteinNames <- readRDS("intersectProteinNames.rds")
  
  result.list <- foreach(i = seq_along(combs.lst)) %dopar% {
    
    vector <- combs.lst[[i]]
    
    statTest <- as.character(vector[1])
    sparcityReduction <- as.character(vector[2])
    normalization <- as.character(vector[3])
    dia <- as.character(vector[4])
    
    print(statTest)
    print(sparcityReduction)
    print(normalization)
    print(dia)
    
    print("---Loading Bootstrap datasets...")
    repList <- readRDS(paste0(dia, ".rds"))
    print("--------------------------")
    subresult.list <- lapply(repList, runAnalysisForEachBootstrapSample, dia=dia, 
                             normalization=normalization, sparcityReduction=sparcityReduction, statTest=statTest,
                             combinedProteinNames=combinedProteinNames, intersectProteinNames=intersectProteinNames)
    
  }
  
  result.list2 <- unlist(result.list, recursive = FALSE)
  result.df <- lapply(result.list2, function(x) x[[2]])
  result.df <- do.call(rbind.data.frame, result.df)
  
  colnames(result.df) <- c("dia", "normalization", "sparcityReduction", "statTest", "groupSize", "nAllProteins", 
                           "nEcoliProteins", "nHumanProteins", "nEcoliProteins.pre", "nHumanProteins.pre", 
                           "medianSampleVariance", "medianProteinVariance", "KS.SignProp", "percNATotal", "percOfRowsWithNAs",
                           "sparcRed.runtime" , "normalization.runtime", "statTest.runtime",
                           "seed.stat",
                           "p.auc", "p.auk", "p.PrecRecallauc",
                           "p.auc.pre", "p.auk.pre", "p.PrecRecallauc.pre", 
                           "p.auc.comb", "p.auk.comb", "p.PrecRecallauc.comb", 
                           "p.auc.intersect", "p.auk.intersect", "p.PrecRecallauc.intersect",
                           "fc.auc", "fc.auc.intersectProtNames",
                           "fc.auc.FCpValRanksum", "fc.auc.FCpValRanksum.intersectProtNames",
                           "RMSEEcoli", "RMSEHuman", "RMSEHumanAndEcoli",
                           "NRMSEEcoli", "NRMSEHuman", "NRMSEHumanAndEcoli",
                           "RMSEEcoli.Intersect", "RMSEHuman.Intersect", "RMSEHumanAndEcoli.Intersect",
                           "NRMSEEcoli.Intersect", "NRMSEHuman.Intersect", "NRMSEHumanAndEcoli.Intersect")
  
  write.csv(result.df, file=paste0("benchmark_results_", batch, "_", batch.size, ".csv"), row.names = TRUE)
  save(result.list, file = paste0("resultlist_", batch,  "_", batch.size, ".Rdata"))
}