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
# SPARSITY REDUCTION

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

# Kolmogorov-Smirnov Test, 
# Returns percentage of significant samples when KS test was applied on single samples compared to all samples combined
kolSmirTestSignProp <- function(mtx) {
  pvals.mtx <- apply(mtx, 2, function(x) stats::ks.test(x, mtx)$p.value)
  cnt <- sum(pvals.mtx<0.05)
  signProportion <- cnt/length(pvals.mtx)
  return(signProportion)
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
  
  percNATotal <- mean(is.na(df)) * 100 
  percOfRowsWithNAs <- sum(apply(df, 1, anyNA))/nrow(df) * 100
  
  sparcRed.runtime <- system.time({ 
    # Sparcity Reduction
    df.sr <- getSparcityReducedDf(sparcityReduction, df) 
  })   
  
  normalization.runtime <- system.time({ 
    # Normalization
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
  
  nEcoli <- nrow(stats.df[which(grepl("ECOLI", stats.df$Protein)),])
  nHuman <- nrow(stats.df[which(grepl("HUMAN", stats.df$Protein)),])
  
  summary <- list(dia, normalization, sparcityReduction, statTest, groupSize=group.size, nAllProteins=nrow(stats.df), 
                  nEcoliProteins=nEcoli, nHumanProteins=nHuman, nEcoliProteins.pre=nEcoli.pre, nHumanProteins.pre=nHuman.pre,
                  medianSampleVariance=medianSampleVariance, medianProteinVariance=medianProteinVariance, KS.SignProp=KS.SignProp, percNATotal=percNATotal, percOfRowsWithNAs=percOfRowsWithNAs,
                  sparcRed.runtime = sparcRed.runtime, normalization.runtime = normalization.runtime, statTest.runtime = statTest.runtime,
                  seed.stat=seed.stat)
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
  statTests <- c("ttest", "limma", "GLMgamma", "SAM", "ROTS") 
  
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
                           "seed.stat")
  
  write.csv(result.df, file=paste0("benchmark_results_", batch, "_", batch.size, ".csv"), row.names = TRUE)
  save(result.list, file = paste0("resultlist_", batch,  "_", batch.size, ".Rdata"))
}