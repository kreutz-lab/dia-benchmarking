library(data.table)
library(rlist)
library(tidyr)
library(dplyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#setwd('..')
#diaWorkflowResults <- lapply("Converted_to_same_format.RData", function(x) mget(load(x)))
diaWorkflowResultsSkyline_DIANNProtSumm <- unlist(lapply("Skyline_DIANNProteinSummarization.RData", function(x) mget(load(x))), recursive = FALSE)
diaWorkflowResultsSkyline_MSstatsProtSumm <- unlist(lapply("Skyline_MSstatsProteinSummarization.RData", function(x) mget(load(x))), recursive = FALSE)

diaWorkflowResultsSkyline_DIANNProtSumm <- list.remove(diaWorkflowResultsSkyline_DIANNProtSumm, ".Random.seed")

names(diaWorkflowResultsSkyline_DIANNProtSumm)[which(grepl("DIANN", names(diaWorkflowResultsSkyline_DIANNProtSumm)))] <-
  sub("DIANN", "DIANN_AI_GPF",names(diaWorkflowResultsSkyline_DIANNProtSumm)[which(grepl("DIANN", names(diaWorkflowResultsSkyline_DIANNProtSumm)))])

names(diaWorkflowResultsSkyline_DIANNProtSumm)[which(grepl("PROSIT", names(diaWorkflowResultsSkyline_DIANNProtSumm)))] <-
  sub("PROSIT", "PROSIT_EDIA",names(diaWorkflowResultsSkyline_DIANNProtSumm)[which(grepl("PROSIT", names(diaWorkflowResultsSkyline_DIANNProtSumm)))])

names(diaWorkflowResultsSkyline_DIANNProtSumm) <- paste0("SkylineDIANN_", sub("_export", "",names(diaWorkflowResultsSkyline_DIANNProtSumm)))

names(diaWorkflowResultsSkyline_MSstatsProtSumm) <- sub("RunLevel", "SkylineMSstats",names(diaWorkflowResultsSkyline_MSstatsProtSumm))




diaWorkflowResults <- unlist(lapply("DIAsoftwareOutputProteinLevel.RData", function(x) mget(load(x))), recursive = FALSE)
diaWorkflowResults <- list.remove(diaWorkflowResults, c(".Random.seed", "Skyline_DIANN_AI_GPF", "Skyline_PROSIT_EDIA_GPF", 
                                                        "Skyline_MaxQuant", "Skyline_MSFragger"))

diaWorkflowResults <- append(diaWorkflowResults, c(diaWorkflowResultsSkyline_DIANNProtSumm, diaWorkflowResultsSkyline_MSstatsProtSumm))
diaWorkflowResults <- diaWorkflowResults[order(names(diaWorkflowResults))]


formatMSTable <- function(x) {
  x <- as.data.frame(x)
  x <- x[, c("Protein.Names", "Run", "PG.Quantity", "Species")]
  colnames(x) <- c("Protein.Names", "Sample", "PG.Quantity", "Species")
  
  if (startsWith(toString(x[1,]$Sample), "LE_")){
    x$Sample <- gsub("LE_", "", x$Sample)
    x$Sample <- sub("_[^_]+$", "", x$Sample)
  }
  x$Sample <- gsub("Lymph_|Lymph_Ecoli_|Ecoli_|Lymphnode_|\\.mzML|mzML", "", x$Sample)
  x$Sample <- gsub("0100", "100", x$Sample)
  
  x <- x[!(x$Species == "synthetic" | x$Species == "usion"),]
  x
}

# UNIFY OUTPUT TABLES
diaWorkflowResults <- lapply(diaWorkflowResults, formatMSTable)

# Change Protein names in OSW_MaxQuant from format Entry to Entry Name
OSW_MaxQuant_unique <- unique(diaWorkflowResults[["OSW_MaxQuant"]]$Protein.Names)

# fwrite(list(OSW_MaxQuant_unique), file = "OSW_MaxQuant_unique.txt")

# Derived translation table from https://www.uniprot.org/uploadlists/
mappedProteins <- read.table("OSW_MaxQuant_unique_mappedProteins_20210204.tsv", sep = "\t", header = TRUE)
diaWorkflowResults[["OSW_MaxQuant"]]$Protein.Names <- mappedProteins$Entry.name[match(diaWorkflowResults[["OSW_MaxQuant"]]$Protein.Names, mappedProteins$Entry)]

# Only keep Entry Name from Protein Name in Skyline_MaxQuant, Skyline_PROSIT_EDIA_GPF, Skyline_DIANN_AI_GPF
for (dia in c("SkylineDIANN_MaxQuant", "SkylineDIANN_PROSIT_EDIA_GPF", "SkylineDIANN_DIANN_AI_GPF",
            "SkylineMSstats_MaxQuant", "SkylineMSstats_PROSIT_EDIA_GPF", "SkylineMSstats_DIANN_AI_GPF",
            "SkylineDIANN_MSFragger", "SkylineMSstats_MSFragger", "Spectronaut_MSFragger")){
  diaWorkflowResults[[dia]]$Protein.Names <- sub(".*\\|", "", diaWorkflowResults[[dia]]$Protein.Names)
  
}

# diaWorkflowResults[["SkylineDIANN_MaxQuant"]]$Protein.Names <- sub(".*\\|", "", diaWorkflowResults[["SkylineDIANN_MaxQuant"]]$Protein.Names)
# diaWorkflowResults[["SkylineDIANN_PROSIT_EDIA_GPF"]]$Protein.Names <- sub(".*\\|", "", diaWorkflowResults[["SkylineDIANN_PROSIT_EDIA_GPF"]]$Protein.Names)
# diaWorkflowResults[["SkylineDIANN_DIANN_AI_GPF"]]$Protein.Names <- sub(".*\\|", "", diaWorkflowResults[["SkylineDIANN_DIANN_AI_GPF"]]$Protein.Names)
# diaWorkflowResults[["SkylineMSstats_MaxQuant"]]$Protein.Names <- sub(".*\\|", "", diaWorkflowResults[["SkylineMSstats_MaxQuant"]]$Protein.Names)
# diaWorkflowResults[["SkylineMSstats_PROSIT_EDIA_GPF"]]$Protein.Names <- sub(".*\\|", "", diaWorkflowResults[["SkylineMSstats_PROSIT_EDIA_GPF"]]$Protein.Names)
# diaWorkflowResults[["SkylineMSstats_DIANN_AI_GPF"]]$Protein.Names <- sub(".*\\|", "", diaWorkflowResults[["SkylineMSstats_DIANN_AI_GPF"]]$Protein.Names)
# 
# diaWorkflowResults[["Spectronaut_MSFragger"]]$Protein.Names <- gsub("'", '', sub(".*\\|", "", diaWorkflowResults[["Spectronaut_MSFragger"]]$Protein.Names))
# 
# diaWorkflowResults[["SkylineDIANN_MSFragger"]]$Protein.Names <- sub(".*\\|", "", diaWorkflowResults[["SkylineDIANN_MSFragger"]]$Protein.Names)
# diaWorkflowResults[["SkylineMSstats_MSFragger"]]$Protein.Names <- sub(".*\\|", "", diaWorkflowResults[["SkylineMSstats_MSFragger"]]$Protein.Names)


addLeadingZeros <- function(el){
  if (length(el) == 1){
    elFinal <- formatC(as.numeric(el),width=3,format='f',digits=0,flag='0')
  } else if (length(el) == 2){
    elZerosAdded <- formatC(as.numeric(el[2]),width=3,format='f',digits=0,flag='0')
    elFinal <- paste0(el[1], "_", elZerosAdded)
  }
  elFinal
}

diaWorkflowResults[["OSW_MaxQuant"]]$Sample <- unlist(lapply( strsplit(unique(diaWorkflowResults[["OSW_MaxQuant"]]$Sample),split='_', fixed=TRUE), addLeadingZeros))
diaWorkflowResults[["OSW_DIANN_AI_GPF"]]$Sample <- unlist(lapply( strsplit(unique(diaWorkflowResults[["OSW_DIANN_AI_GPF"]]$Sample),split='_', fixed=TRUE), addLeadingZeros))
diaWorkflowResults[["OSW_MSFragger"]]$Sample <- unlist(lapply( strsplit(unique(diaWorkflowResults[["OSW_MSFragger"]]$Sample),split='_', fixed=TRUE), addLeadingZeros))


generateWideFormattedDf <- function(df) {
  # remove duplicates from multiple precursors carrying the same quantitative information for a certain protein
  df <- dplyr::distinct(df)

  df <- reshape(df, idvar="Protein.Names", timevar="Sample",  direction="wide")
  colnames(df) <- sub("PG.Quantity.", "", colnames(df))

  row.names(df) <- df$Protein.Names
  df$Protein.Names <- NULL
  df[df == 0] <- NA
  df <- log2(df)

  #df <- df[,grepl("1-6", colnames(df)) | grepl("1-12", colnames(df)) | grepl("1-25", colnames(df))]
  
  # For "Skyline_DIANN_AI_GPF", "Skyline_MaxQuant",  "Skyline_PROSIT_EDIA_GPF" "1-6_028" is missing, 
  # add empty column with that name as a placeholder
  if (!"1-6_028" %in% names(df)){
    df["1-6_028"] <- NA
  }
  
  return(df)
}


cleanData <- function(x) {
  # remove proteins with both human and e.coli annotation
  x <- x[!(grepl("HUMAN", x$Protein.Names) & grepl("ECOLI", x$Protein.Names)),]
  
  x <- x[with(x, order(Protein.Names, Sample)),]
  x$Species <- NULL
  # If there are duplicates in terms of Protein.Names and Samples combined, check if one them is NA for PG.Quantity and remove respective entry
  if (sum(duplicated(x[,c('Protein.Names', 'Sample')]))>0){
    frequencies <- x %>%
      group_by(Protein.Names, Sample) %>%
      dplyr::summarise(Freq = n())
    
    x.freqAdded <- merge(x, frequencies)
    noDuplicate <- x.freqAdded$Freq == 1
    missing <- is.na(x$PG.Quantity)
    x <- x[noDuplicate | (!noDuplicate & !missing), ]
  }
  x <- generateWideFormattedDf(x)
  
  # Remove empty rows
  x <- x[rowSums(is.na(x)) != ncol(x),]
  x
}

diaWorkflowResults <- lapply(diaWorkflowResults, cleanData)

# cols <- colnames(diaWorkflowResults[[1]])
# design <- data.frame(colnames=cols, dilution=cols)
# 
# #design[!(grepl("1-6", design$colnames) | grepl("1-12", design$colnames) | grepl("1-25", design$colnames)), ]$dilution <- "undiluted"
# design[grepl("1-6", design$colnames), ]$dilution <- "X6"
# design[grepl("1-12", design$colnames), ]$dilution <- "X12"
# design[grepl("1-25", design$colnames), ]$dilution <- "X25"
# row.names(design) <- design$colnames
# 
# design$dilution <- factor(design$dilution,
#                           levels = c("X25", "X12", "X6"))
# design <- design[order(design$dilution), ]
# design$colnames <- NULL


getDesign <- function(df) {
  cols <- colnames(df)
  design <- data.frame(colnames=cols, dilution=cols)
  
  design[!(grepl("1-6", design$colnames) | grepl("1-12", design$colnames) | grepl("1-25", design$colnames)), ]$dilution <- "Lymph nodes"
  design[grepl("1-6", design$colnames), ]$dilution <- "Lymph nodes + 1:6 E.coli"
  design[grepl("1-12", design$colnames), ]$dilution <- "Lymph nodes + 1:12 E.coli"
  design[grepl("1-25", design$colnames), ]$dilution <- "Lymph nodes + 1:25 E.coli"
  row.names(design) <- design$colnames
  
  design$dilution <- factor(design$dilution,
                            levels = c("Lymph nodes", "Lymph nodes + 1:25 E.coli", "Lymph nodes + 1:12 E.coli", "Lymph nodes + 1:6 E.coli"))
  
  design <- design[order(design$dilution), ]
  design$colnames <- NULL
  design
}

design <- getDesign(diaWorkflowResults[[1]])

samplenames <- row.names(design)
idx25 <- which(design$dilution == "Lymph nodes + 1:25 E.coli")
idx12 <- which(design$dilution == "Lymph nodes + 1:12 E.coli")

saveRDS(design, file = "design.rds")

diaWorkflowResults <- lapply(diaWorkflowResults,function(x, samplenames, idx25, idx12) {
  # Bring to same order as design data frame
  x <- x[, samplenames]
  x
}, samplenames=samplenames)

saveRDS(diaWorkflowResults, file = "diaWorkflowResults_allDilutions.rds")


diaWorkflowResults <- lapply(diaWorkflowResults,function(x, samplenames, idx25, idx12) {

  # Keep only samples of dilution 12 or 25
  x <- x[, c(idx25, idx12)]

  # Remove empty rows
  x <- x[rowSums(is.na(x)) != ncol(x),]

  samplenames <- samplenames[grepl("1-25", samplenames) | grepl("1-12", samplenames)]
  x <- x[,samplenames]
  x
}, samplenames=samplenames, idx25=idx25, idx12 = idx12)

design$dilution2 <- design$dilution
design_1225only <- design[grepl("1:25", design$dilution) | grepl("1:12", design$dilution),]
idx25_1225only <- which(design_1225only$dilution == "Lymph nodes + 1:25 E.coli")
idx12_1225only <- which(design_1225only$dilution == "Lymph nodes + 1:12 E.coli")


set.seed(50)
indices <- list()
# 100 iterations
for (i in seq(1:100)){
  # Sample sizes from 3 to 23
  for(sampleSize in seq(3, length(idx25_1225only))){
    indices25 <- sample(idx25_1225only, sampleSize, replace = T)
    indices12 <- sample(idx12_1225only, sampleSize, replace = T)
    indices <- list.append(indices, c(indices25, indices12))
  }
}

# # same sample indices are used for all DIA workflows
# repList <- list()
# repList <- lapply(diaWorkflowResults,function(x, repList, indices) {
#   lst <- lapply(indices,function(idx, x) {
#     x <- x[,idx]
#     x
#   },  x=x)
# 
#   repList <- list.append(repList, lst)
#   repList
# }, repList=repList, indices=indices)
# 
# 
# repList <- unlist(repList, recursive = FALSE)


# Benchmark analysis is only conducted on Skyline data with protein summarization via MSstats and not via DIANN
diaWorkflowResults <- diaWorkflowResults[!names(diaWorkflowResults) %in% c("SkylineDIANN_DIANN_AI_GPF", "SkylineDIANN_MaxQuant", 
                                                                           "SkylineDIANN_MSFragger", "SkylineDIANN_PROSIT_EDIA_GPF")]
names(diaWorkflowResults) <- sub("SkylineMSstats", "Skyline",names(diaWorkflowResults))

combinedProteinNames <- c()
combinedProteinNames <- lapply(diaWorkflowResults, function(x, combinedProteinNames){
  names <- unique(unlist(strsplit(row.names(x), ";")))
  combinedProteinNames <- c(combinedProteinNames, names)
  combinedProteinNames
}, combinedProteinNames=combinedProteinNames)

# proteins appearing in at least 14/17 (82%) of the DIA workflows (= 'core dataset')
intersectProteinNames <- data.frame(table(unlist(combinedProteinNames)))
intersectProteinNames <- intersectProteinNames[intersectProteinNames$Freq>=14,]$Var1

# proteins appearing in at least one DIA workflow
combinedProteinNames <- sort(unique(unlist(combinedProteinNames)))

#save(diaWorkflowResults, indices, intersectProteinNames, combinedProteinNames, file="Converted_to_same_format_wide_plusIndices_plusIntersectingAndCombinedProteinNames.RData")
# save(diaWorkflowResults, indices, intersectProteinNames, combinedProteinNames, file="DIAsoftwareOutputProteinLevel_1to12And1to25Only_wideFormat_withBootstrapIndicesAndIntersectAndCombinedProteinNames.RData")

saveRDS(intersectProteinNames, file = "intersectProteinNames.rds")
saveRDS(combinedProteinNames, file = "combinedProteinNames.rds")
saveRDS(indices, file = "indices.rds")
saveRDS(diaWorkflowResults, file = "diaWorkflowResults.rds")

session <- sessionInfo()
sink(paste0("sessionInfo_bootstrapping.txt"))
print(session)
sink()
