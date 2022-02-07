library(dplyr)
library(MBQN)
library(matrixTests)
library(limma)
library(qvalue)
library(genefilter)
library(samr)
# BiocManager::install("ROTS")
library(ROTS)
library(parallel)
library(DescTools)
library(matrixcalc)
library(psych) 
library(doParallel)

library(stats)
library(pcaMethods)
library(foreach)
library(rlist)
library(matrixStats) 

#################################################################################
# NORMALIZATION
# "unnormalized", "TRQN", "QN", "median"
getNormalizedDf <- function(modus, df) {
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
# SPARSITY REDUCTION

getSparsityReducedDf <- function(modus, df) {
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
    if (modelType == "glm") {
      # res <- as.vector(coef(summary(glm(x2 ~ y, family  = Gamma(link = "identity"), data = tmp)))[,"Pr(>|t|)"][2])
      res <- as.vector(coef(summary(glm(x2 ~ y, family  = Gamma(link = "log"), data = tmp)))[,"Pr(>|t|)"][2]) # (link = "identity") was wrongly added
      #res <- as.vector(coef(summary(glm(x2 ~ y, family  = "Gamma", data = tmp)))[,"Pr(>|t|)"][2]) # (link = "identity") was wrongly added
      
    }
    # if (modelType == "lm") {
    #   res <- as.vector(coef(summary(lm(x2 ~ y, data=tmp)))[,"Pr(>|t|)"][2])
    # } else if (modelType == "glm") {
    #   res <- as.vector(coef(summary(glm(x2 ~ y, family  = Gamma(link = "identity"), data = tmp)))[,"Pr(>|t|)"][2])
    # } else if (modelType == "lasso") {
    #   res <- as.vector(coef(summary(lasso2::l1ce(x2 ~ y, data=tmp)))[,"Pr(>|Z|)"][2])
    # }
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
  if (modus == "ttestVarEqual") {
    output <- matrixTests::row_t_equalvar(df[, 1:group.size], df[, (group.size+1):(2*group.size)], alternative = "two.sided", mu = 0, conf.level = 0.95)
    stats.df <- data.frame(Protein = row.names(df), pValue = output$pvalue)
  } else if (modus == "ttestVarNotEqual") {
    output <- matrixTests::row_t_welch(df[, 1:group.size], df[, (group.size+1):(2*group.size)], alternative = "two.sided", mu = 0, conf.level = 0.95)
    stats.df <- data.frame(Protein = row.names(df), pValue = output$pvalue)
  } else if (modus == "Wilcoxon") {
    output <- matrixTests::row_wilcoxon_twosample(df[, 1:group.size], df[, (group.size+1):(2*group.size)], alternative = "two.sided", mu = 0, exact = NA, correct = TRUE)
    stats.df <- data.frame(Protein = row.names(df), pValue = output$pvalue)
  # } else if (modus == "RankProduct") {
  #   seed <- sample(1:1000000000, 1)
  #   groups <-  c(rep(0, group.size), rep(1, group.size))
  #   output <- RankProd::RankProducts(df, groups, rand=seed)
  #   stats.df <- data.frame(Protein = row.names(df), pValue = apply(output$pval,1,min))
  } else if (modus == "limma") {
    # limma
    # library(limma)
    design <- model.matrix(~groups)
    fit <- lmFit(df, design)
    fit2 <- eBayes(fit)
    output <- topTable(fit2, coef = 2, number = nrow(df), sort.by = "none", adjust.method="BH")
    stats.df <- data.frame(Protein = row.names(df), pValue = output$P.Value)
  # } else if (modus == "ttest") {
  #   pVals <- rowttests(data.matrix(df),as.factor(groups), na.rm = T)$p.value
  #   stats.df <- data.frame(Protein = row.names(df), pValue = pVals)
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
  } else if (modus == "GLMgamma") {
    # GLM-Gamma
    seed <- sample(1:1000000000, 1)
    set.seed(seed)
    exp.df <- 2^df # Gamma regression needs to be done on non-lol transformed data
    stats.lst <- apply(exp.df, 1, multiModel, group.size=group.size, modelType="glm")
    stats.df <- statsLstToDf(stats.lst)
  } else {
    print("Undefined modus")
  }
  return(list(stats.df, seed))
}

#################################################################################
# DATA CHARACTERISTICS


# Kolmogorov-Smirnov Test, 
# Returns percentage of significant samples when KS test was applied on single samples compared to all samples combined
kolSmirTestSignProp <- function(mtx) {
  pvals.mtx <- apply(mtx, 2, function(x) stats::ks.test(x, mtx)$p.value)
  cnt <- sum(pvals.mtx<0.05)
  signProportion <- cnt/length(pvals.mtx)
  return(signProportion)
}

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

getMoreCharacteristics <- function(mtx, withNAs=TRUE){
   KS.SignProp <- kolSmirTestSignProp(mtx)
  
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
    resultvec <- c(KS.SignProp = KS.SignProp,
                   entropy = entropy,
                   kurtosis = kurtosis, 
                   meanDeviation = meanDeviation,
                   skewness = skewness,
                   uniformity = uniformity,
                   variance = variance, 
                   RMS = RMS,
                   var.groups.ratio = var.groups.ratio)
  } else {
    t.mtx <- t(mtx)
    t.mtx <- t.mtx[ , which(apply(t.mtx, 2, var) != 0)] # Remove zero variance columns 
    pca <- stats::prcomp(t.mtx, scale.=T)
    eigs <- pca$sdev^2
    prctPC1 <- eigs[1]/sum(eigs)
    prctPC2 <- eigs[2]/sum(eigs)
    
    elongation <- sqrt(eigs[2] / eigs[1]) # elongation
    flatness <- sqrt(eigs[length(eigs)]/eigs[1]) # flatness
    resultvec <- c(KS.SignProp = KS.SignProp,
                   entropy = entropy,
                   kurtosis = kurtosis, 
                   meanDeviation = meanDeviation,
                   skewness = skewness,
                   uniformity = uniformity,
                   variance = variance, 
                   RMS = RMS,
                   var.groups.ratio = var.groups.ratio,
                   prctPC1 = prctPC1, 
                   prctPC2 = prctPC2,
                   elongation = elongation,
                   flatness = flatness)
  }
  return(resultvec)
}

runBPTest <- function(x, group.size){
  tmp <- data.frame(x2=unlist(as.vector(x)), y=as.factor(c(rep("X25", group.size), rep("X12", group.size))))
  m <- stats::lm(x2 ~ y, data=tmp)
  bp.res <- lmtest::bptest(m, studentize = TRUE) 
  bp.res[["p.value"]]
}


getDataCharacteristics <- function(df) {
  mtx <- as.matrix(df)
  medianSampleVariance <- median(apply(df, 2, var,  na.rm=TRUE))
  medianProteinVariance <- median(unname(apply(df, 1, var,  na.rm=TRUE)), na.rm = TRUE)
  #KS.SignProp <- kolSmirTestSignProp(as.matrix(df))
  percNATotal <- mean(is.na(df)) * 100 
  percOfRowsWithNAs <- sum(apply(df, 1, anyNA))/nrow(df) * 100
  
  #percNATotal2 <- mean(is.na(mtx)) * 100 
  
  characts.wNAs <- getMoreCharacteristics(mtx, withNAs=TRUE)
  # names(characts.wNAs) <- c("entropy.wNAs",
  #                           "kurtosis.wNAs", 
  #                           "meanDeviation.wNAs",
  #                           "skewness.wNAs",
  #                           "uniformity.wNAs",
  #                           "variance.wNAs", 
  #                           "RMS.wNAs",
  #                           "var.groups.ratio.wNAs")
  names(characts.wNAs) <- paste0(names(characts.wNAs), ".wNAs")
  
  # "prctPC1.wNAs", 
  #"elongation.wNAs",
  #"flatness.wNAs")
  
  mtx <- mtx[rowSums(is.na(mtx)) == 0, ]
  nProteins.woNAs <- nrow(mtx) # number of proteins with no NAs
  characts.woNAs <- getMoreCharacteristics(mtx, withNAs=FALSE)
  # names(characts.woNAs) <- c("entropy.woNAs",
  #                            "kurtosis.woNAs", 
  #                            "meanDeviation.woNAs",
  #                            "skewness.woNAs",
  #                            "uniformity.woNAs",
  #                            "variance.woNAs", 
  #                            "RMS.woNAs",
  #                            "var.groups.ratio.woNAs",
  #                            "prctPC1.woNAs", 
  #                            "elongation.woNAs",
  #                            "flatness.woNAs")
  names(characts.woNAs) <- paste0(names(characts.woNAs), ".woNAs")
  
  group.size <- ncol(mtx)/2
  BPTest.lst <- apply(mtx, 1, runBPTest, group.size=group.size)
  
  # heterosc.woNAs <- sum(BPTest.lst < 0.05) / length(BPTest.lst)
  
  qobj <- qvalue::qvalue(p = BPTest.lst)
  heterosc.oneMinuspi0 <- 1 - qobj$pi0
  
  datacharacts <- c(medianSampleVariance = medianSampleVariance, medianProteinVariance = medianProteinVariance, 
                    percNATotal = percNATotal, percOfRowsWithNAs = percOfRowsWithNAs, 
                    characts.wNAs, characts.woNAs, heterosc.oneMinuspi0=heterosc.oneMinuspi0,
                    nProteins.woNAs=nProteins.woNAs)
}

getNumberOfProteins <- function(nEcoli.pre, nHuman.pre, stats.df, combinedProteinNames, intersectProteinNames, DiaWorkflowProteinNames) {
  
  nEcoli <- nrow(stats.df[which(grepl("ECOLI", stats.df$Protein)),])
  nHuman <- nrow(stats.df[which(grepl("HUMAN", stats.df$Protein)),])
  
  nEcoli.comb <- length(combinedProteinNames[grepl("_ECOLI", combinedProteinNames)])
  nHuman.comb <- length(combinedProteinNames[grepl("_HUMAN", combinedProteinNames)])

  nEcoli.intersect <- length(intersectProteinNames[grepl("_ECOLI", intersectProteinNames)])
  nHuman.intersect <- length(intersectProteinNames[grepl("_HUMAN", intersectProteinNames)])

  nEcoli.diaWorkflow <- length(DiaWorkflowProteinNames[grepl("_ECOLI", DiaWorkflowProteinNames)])
  nHuman.diaWorkflow <- length(DiaWorkflowProteinNames[grepl("_HUMAN", DiaWorkflowProteinNames)])
  
  c(nEcoli.pre=nEcoli.pre, nHuman.pre=nHuman.pre,
    nEcoli=nEcoli, nHuman=nHuman, 
    nEcoli.comb = nEcoli.comb, nHuman.comb = nHuman.comb,
    nEcoli.intersect = nEcoli.intersect, nHuman.intersect = nHuman.intersect,
    nEcoli.diaWorkflow = nEcoli.diaWorkflow, nHuman.diaWorkflow = nHuman.diaWorkflow)
}

#################################################################################
# pAUC calculation
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
    } else {
      diffs.x <- sp[-1] - sp[-length(sp)]
      means.vert <- (se[-1] + se[-length(se)])/2
      auc <- sum(means.vert * diffs.x)
    }
  } else { # Partial area
    if (partial.auc.focus == "sensitivity") {
      # if we focus on SE, just swap and invert x and y and the computations for SP will work
      x <- rev(se)
      y <- rev(sp)
    }else {
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
    } else if (percent) {
      auc <- (100+((auc-min)*100/(max-min)))/2 # McClish formula adapted for %
    } else {
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

getpValCurvedf <- function(stats.df, nEcoli.pre=NA, nHuman.pre=NA) {
  
  if (is.na(nEcoli.pre)){
    totalEcoli <- nrow(stats.df[which(grepl("ECOLI", row.names(stats.df))),])
    totalHuman <- nrow(stats.df[which(grepl("HUMAN", row.names(stats.df))),])
  } else {
    totalEcoli <- nEcoli.pre
    totalHuman <- nHuman.pre
  }
  
  #p.curve.lst <- lapply(unique(c(0, sort(stats.df$pValue)[c(TRUE, rep(FALSE, times=9))], 1)), function(pvalue, stats.df){
  p.curve.lst <- lapply(unique(c(0, sort(stats.df$pValue), 1)), function(pvalue, stats.df){
    
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
  partial.auc.corrects <- c(FALSE) # c(TRUE, FALSE)
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
    if (length(auc) == 0) auc <- NA # length 0 can happen if 1-specificity for pAUC can't be reached, e.g due too few proteins being left after sparsity reduction
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


getStatsProteinNames <- function(combIntersect, proteinNames, stats.df, dia=NULL) {
  # if (combIntersect =="combined"){
  #   proteinNames <- readRDS("combinedProteinNames.rds")
  # } else if (combIntersect =="intersect"){
  #   proteinNames <- readRDS("intersectProteinNames.rds")
  # } else if (combIntersect =="diaWorkflow"){
  #   proteinNames <- readRDS(paste0(dia, "_ProteinNames.rds"))
  # }
  nEcoli <- length(proteinNames[grepl("_ECOLI", proteinNames)])
  nHuman <- length(proteinNames[grepl("_HUMAN", proteinNames)])
  
  if (combIntersect %in% c("combined", "intersect")){
    intersectBool <- apply(stats.df, 1, function(x) {
      length(intersect(unlist(base::strsplit(x[1], ";")), proteinNames))>0
    })
    stats.df.protNames <- stats.df[intersectBool, ]
  }  else if (combIntersect =="diaWorkflow"){
    stats.df.protNames <- stats.df
  }
  
  list(stats.df.protNames=stats.df.protNames, nEcoli=nEcoli, nHuman=nHuman)
}

getRegValsPropAndPAucs <- function(combIntersect=c("combined", "intersect", "diaWorkflow"), proteinNames, stats.df, dia) {
  row.names(stats.df) <- stats.df$Protein
  stats.protNames <- getStatsProteinNames(combIntersect = combIntersect, proteinNames=proteinNames, stats.df, dia)
  regpValsProp <- getRegulatedProp(pValueVec=stats.protNames[["stats.df.protNames"]]$pValue) 
  stats.protNames[["stats.df.protNames"]]$pValue[is.na(stats.protNames[["stats.df.protNames"]]$pValue)] <- 1 # Replace NAs with pValue of 1
  
  paucs <- getPartialAUCsResults(stats.df = stats.protNames[["stats.df.protNames"]], 
                                 nEcoli = stats.protNames[["nEcoli"]], 
                                 nHuman = stats.protNames[["nHuman"]])
  
  res <- c(paucs, 
           regpValsProp=regpValsProp)
  
  names(res) <- paste0(names(res), ".", combIntersect)
  
  return(res)
}

#################################################################################

# RMSE

rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2,na.rm = TRUE))
}

getRMSE <- function(stats.df) {
  stats.dfEcoli <- stats.df[grepl("ECOLI", stats.df$Protein), ]
  stats.dfHuman <- stats.df[grepl("HUMAN", stats.df$Protein), ]
  
  Ecoli <- rmse(actual=stats.dfEcoli$log2FC, predicted=stats.dfEcoli$log2FCPredicted)
  Human <- rmse(actual=stats.dfHuman$log2FC, predicted=stats.dfHuman$log2FCPredicted)
  HumanAndEcoli <- rmse(actual=stats.df$log2FC, predicted=stats.df$log2FCPredicted)
  
  list(Ecoli, Human, HumanAndEcoli)
}

getRMSEResults <- function(combIntersect=c("intersect", "diaWorkflow"), proteinNames, stats.df, dia) {
  
  stats.protNames <- getStatsProteinNames(combIntersect = combIntersect, proteinNames = proteinNames, stats.df, dia)
  
  stats.df2 <- stats.protNames[["stats.df.protNames"]]
  RMSEs <- getRMSE(stats.df2)
  RMSEEcoli <- RMSEs[[1]]
  RMSEHuman <- RMSEs[[2]]
  RMSEHumanAndEcoli <- RMSEs[[3]]
  
  res <- c(RMSEEcoli=RMSEEcoli, 
           RMSEHuman=RMSEHuman, 
           RMSEHumanAndEcoli=RMSEHumanAndEcoli)
  
  names(res) <- paste0(names(res), ".", combIntersect)
  return(res)
}

#################################################################################

runAnalysisForEachBootstrapSample <- function(bootstrap.dataset, indices, diaWorkflowResults.selected, dia, normalization, sparsityReduction, statTest,
                                              combinedProteinNames, intersectProteinNames, DiaWorkflowProteinNames) {
  df <- diaWorkflowResults.selected[, unlist(indices[bootstrap.dataset])]
  #df <- repList[[bootstrap.dataset]]
  
  # Remove empty rows
  df <- df[rowSums(is.na(df)) != ncol(df), ]     
  
  nEcoli.pre <- sum(grepl("ECOLI", row.names(df)))
  nHuman.pre <- sum(grepl("HUMAN", row.names(df)))  
  group.size <- ncol(df)/2

  data.characts <- getDataCharacteristics(df)
  
  sparsRed.runtime <- system.time({ 
    # Sparsity Reduction
    df.sr <- getSparsityReducedDf(sparsityReduction, df) 
  })   
  
  normalization.runtime <- system.time({ 
    # Normalization
    df <- getNormalizedDf(normalization, df.sr)
  })
  
  log2FC.df <- data.frame(Protein=row.names(df), log2FC=rowMeans(df[, (group.size+1):ncol(df)], na.rm = TRUE)-rowMeans(df[, 1:group.size], na.rm = TRUE))
  
  statTest.runtime <- system.time({ 
    # print("Stats")
    stats.df.lst <- getStatTestResultDf(statTest, df) 
    stats.df <- stats.df.lst[[1]]
    seed.stat <- stats.df.lst[[2]]
  })
  
  sparsRed.runtime <- sparsRed.runtime[["user.self"]]
  normalization.runtime <- normalization.runtime[["user.self"]]
  statTest.runtime <- statTest.runtime[["user.self"]]
  
  stats.df <- merge(stats.df, log2FC.df, by=c("Protein"), all.x=FALSE)
  stats.df$log2FCPredicted <- NA
  stats.df[grepl("ECOLI", stats.df$Protein), ]$log2FCPredicted <- log2(0.24038462/0.11076923) #  1.11778738, was previously wrongly assumed to be log2((1/12)/(1/25))=1.058894 
  stats.df[grepl("HUMAN", stats.df$Protein), ]$log2FCPredicted <- log2(1)
  
  if (sum(is.nan(stats.df$pValue)) > 0) stats.df[is.nan(stats.df$pValue),]$pValue <- NA
  if (sum(is.nan(stats.df$log2FC)) > 0) stats.df[is.nan(stats.df$log2FC),]$log2FC <- NA
  
  numberOfProteins <- getNumberOfProteins(nEcoli.pre, nHuman.pre, stats.df, combinedProteinNames,
                                          intersectProteinNames, DiaWorkflowProteinNames)
  
  
  regValsPropAndPAucs.comb <- getRegValsPropAndPAucs(combIntersect = "combined", proteinNames = combinedProteinNames, stats.df = stats.df, dia=NULL)
  regValsPropAndPAucs.intersect <- getRegValsPropAndPAucs(combIntersect = "intersect", proteinNames = intersectProteinNames, stats.df = stats.df, dia=NULL)
  regValsPropAndPAucs.diaWorkflow <- getRegValsPropAndPAucs(combIntersect = "diaWorkflow", proteinNames = DiaWorkflowProteinNames, stats.df = stats.df, dia=dia)

  RMSE.intersect <- getRMSEResults(combIntersect="intersect", proteinNames = intersectProteinNames, stats.df, dia=NULL) 
  RMSE.diaWorkflow <- getRMSEResults(combIntersect="diaWorkflow", proteinNames = DiaWorkflowProteinNames, stats.df, dia=dia)
  
  summary <- rlist::list.flatten(list(bootstrap.dataset = bootstrap.dataset, dia = dia, normalization = normalization, sparsityReduction = sparsityReduction, statTest = statTest, groupSize=group.size, 
                  as.list(data.characts), 
                  nAllProteins=nrow(stats.df), 
                  as.list(numberOfProteins),
                  sparsRed.runtime = sparsRed.runtime, normalization.runtime = normalization.runtime, statTest.runtime = statTest.runtime,
                  seed.stat=seed.stat,
                  as.list(regValsPropAndPAucs.comb), as.list(regValsPropAndPAucs.intersect), as.list(regValsPropAndPAucs.diaWorkflow),
                  as.list(RMSE.intersect), as.list(RMSE.diaWorkflow)))
  return(list(stats.df, summary))
}

################################################################################

# batch <- 1 # max 32 bei batch.size 20, max 64 bei batch.size 10, max 128 bei batch.size 5
run <- function(batch, batch.size) {
  
  # GET ALL PARAMETER COMBINATIONS
  dias <- sort(c("DIANN_DIANN_AI", "DIANN_DIANN_AI_GPF", "DIANN_MaxQuant", "DIANN_MSFragger", 
                 "DIANN_PROSIT_EDIA_GPF", "OSW_DIANN_AI_GPF", "OSW_MaxQuant", 
                 "OSW_MSFragger", "Skyline_DIANN_AI_GPF", 
                 "Skyline_MaxQuant", "Skyline_MSFragger", "Skyline_PROSIT_EDIA_GPF", 
                 "Spectronaut_DIANN_AI_GPF", "Spectronaut_DirectDIA", "Spectronaut_MaxQuant", 
                 "Spectronaut_MSFragger", "Spectronaut_PROSIT_EDIA_GPF"))
  normalizations <- c("unnormalized", "TRQN", "QN", "median")
  sparsityReductions <- c("NoSR", "SR66", "SR90")
  statTests <- c("ttestVarEqual", "ttestVarNotEqual", "GLMgamma", "limma", "Wilcoxon", "SAM", "ROTS") 
  
  combs <- expand.grid(statTests, sparsityReductions, normalizations, dias)
  colnames(combs) <- c("statTest", "sparsityReduction", "normalization", "dia")
  
  if (batch*batch.size > nrow(combs)){
    combs <- combs[(((batch-1)*batch.size)+1):nrow(combs),]
  } else{
    combs <- combs[(((batch-1)*batch.size)+1):(batch*batch.size),]
  }
  
  combs.lst <- as.list(as.data.frame(t(combs)))
  

  
  combinedProteinNames <- readRDS("combinedProteinNames.rds")
  intersectProteinNames <- readRDS("intersectProteinNames.rds")
  # load("/Users/admin/Desktop/PhD/202110_dia-benchmarking_rerun/DIAsoftwareOutputProteinLevel_1to12And1to25Only_wideFormat_withBootstrapIndicesAndIntersectAndCombinedProteinNames.RData")
  # saveRDS(diaWorkflowResults, file = "diaWorkflowResults.rds")
  indices <- readRDS("indices.rds")
  diaWorkflowResults <- readRDS("diaWorkflowResults.rds")
  
  #result.list <- foreach(i = seq_along(combs.lst)) %dopar% {
  #result.list <- foreach(i = seq_along(combs.lst)) %do% {  
    vector <- combs.lst[[1]]
    
    statTest <- as.character(vector[1])
    sparsityReduction <- as.character(vector[2])
    normalization <- as.character(vector[3])
    dia <- as.character(vector[4])
    
    print(statTest)
    print(sparsityReduction)
    print(normalization)
    print(dia)
    
    #DiaWorkflowProteinNames <- readRDS(paste0(dia, "_ProteinNames.rds"))
    
    # print("---Loading Bootstrap datasets...")
    # repList <- readRDS(paste0(dia, ".rds"))
    print("--------------------------")
    # subresult.list <- lapply(seq_along(repList), runAnalysisForEachBootstrapSample, repList=repList, dia=dia, 
    #                          normalization=normalization, sparsityReduction=sparsityReduction, statTest=statTest,
    #                          combinedProteinNames=combinedProteinNames, intersectProteinNames=intersectProteinNames, DiaWorkflowProteinNames=DiaWorkflowProteinNames)
    diaWorkflowResults.selected <- diaWorkflowResults[[dia]]
    DiaWorkflowProteinNames <- row.names(diaWorkflowResults.selected)
    
    # result.list <- lapply(seq_along(indices), runAnalysisForEachBootstrapSample, indices=indices, diaWorkflowResults.selected=diaWorkflowResults.selected, dia=dia, 
    #                          normalization=normalization, sparsityReduction=sparsityReduction, statTest=statTest,
    #                          combinedProteinNames=combinedProteinNames, intersectProteinNames=intersectProteinNames, DiaWorkflowProteinNames=DiaWorkflowProteinNames)
    # 
    
    procs <- as.numeric(Sys.getenv("SLURM_NTASKS"))
    registerDoParallel(cores=procs)
    
    result.list <- foreach(i = seq_along(indices)) %dopar% {
      runAnalysisForEachBootstrapSample(i, indices=indices, diaWorkflowResults.selected=diaWorkflowResults.selected, dia=dia, 
      normalization=normalization, sparsityReduction=sparsityReduction, statTest=statTest,
      combinedProteinNames=combinedProteinNames, intersectProteinNames=intersectProteinNames, DiaWorkflowProteinNames=DiaWorkflowProteinNames)
    }
  #}
  
  #result.list2 <- unlist(result.list, recursive = FALSE)
  result.df <- lapply(result.list, function(x) x[[2]])
  
  result.df <- do.call(rbind.data.frame, result.df)

  # colnames(result.df) <- c("bootstrap.dataset", "dia", "normalization", "sparsityReduction", "statTest", "groupSize", "nAllProteins", 
  #                          "nEcoliProteins", "nHumanProteins", "nEcoliProteins.pre", "nHumanProteins.pre", 
  #                          "medianSampleVariance", "medianProteinVariance", "KS.SignProp", "percNATotal", "percOfRowsWithNAs",
  #                          "sparsRed.runtime" , "normalization.runtime", "statTest.runtime",
  #                          "seed.stat")
  #colnames(result.df) <- names(subresult.list[[1]][[2]])
 
  if (batch.size == 1){
    session <- sessionInfo()
    sink(paste0("sessionInfo_", batch, "_",batch.size, 
                "_", dia, "_", normalization, "_", sparsityReduction, "_", statTest,".txt"))
    print(session)
    sink()
    
    write.csv(result.df, file=paste0("benchmark_results_", batch, "_", batch.size, 
                                     "_", dia, "_", normalization, "_", sparsityReduction, "_", statTest,".csv"), row.names = FALSE)
    save(result.list, file = paste0("resultlist_", batch,  "_", batch.size, 
                                    "_", dia, "_", normalization, "_", sparsityReduction, "_", statTest, ".Rdata"))
  } else {
    session <- sessionInfo()
    sink(paste0("sessionInfo_", batch, "_",batch.size, ".txt"))
    print(session)
    sink()
    
    write.csv(result.df, file=paste0("benchmark_results_", batch, "_", batch.size, ".csv"), row.names = FALSE)
    save(result.list, file = paste0("resultlist_", batch,  "_", batch.size, ".Rdata"))
    
  }
  
}

