library(data.table)
library(rlist)
library(data.table) # install if not installed already
library(tidyr)
library(dplyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')
#diaWorkflowResults <- lapply("Converted_to_same_format.RData", function(x) mget(load(x)))
diaWorkflowResults <- lapply("DIAsoftwareOutputProteinLevel.RData", function(x) mget(load(x)))

diaWorkflowResults <- unlist(diaWorkflowResults,recursive=FALSE)

diaWorkflowResults <- list.remove(diaWorkflowResults, ".Random.seed")

# UNIFY OUTPUT TABLES
diaWorkflowResults <- lapply(diaWorkflowResults, function(x) {
  x <- as.data.frame(x)
  x <- x[, c("Protein.Names", "Run", "PG.Quantity", "Species")]
  colnames(x) <- c("Protein.Names", "Sample", "PG.Quantity", "Species")

  if (startsWith(toString(x[1,]$Sample), "LE_")){
    x$Sample <- gsub("LE_", "", x$Sample)
    x$Sample <- sub("_[^_]+$", "", x$Sample)
  }
  x$Sample <- gsub("Lymph_|Lymph_Ecoli_|Ecoli_|Lymphnode_|.mzML", "", x$Sample)
  x$Sample <- gsub("0100", "100", x$Sample)

  x <- x[!(x$Species=="synthetic"),]
  x
})

# Change Protein names in OSW_MaxQuant from format Entry to Entry Name
OSW_MaxQuant_unique <- unique(diaWorkflowResults[["OSW_MaxQuant"]]$Protein.Names)

# fwrite(list(OSW_MaxQuant_unique), file = "OSW_MaxQuant_unique.txt")

# Derived translation table from https://www.uniprot.org/uploadlists/
mappedProteins <- read.table("OSW_MaxQuant_unique_mappedProteins_20210204.tsv", sep = "\t", header = TRUE)
diaWorkflowResults[["OSW_MaxQuant"]]$Protein.Names <- mappedProteins$Entry.name[match(diaWorkflowResults[["OSW_MaxQuant"]]$Protein.Names, mappedProteins$Entry)]

# Only keep Entry Name from Protein Name in Skyline_MaxQuant, Skyline_PROSIT_EDIA_GPF, Skyline_DIANN_AI_GPF
diaWorkflowResults[["Skyline_MaxQuant"]]$Protein.Names <- sub(".*\\|", "", diaWorkflowResults[["Skyline_MaxQuant"]]$Protein.Names)
diaWorkflowResults[["Skyline_PROSIT_EDIA_GPF"]]$Protein.Names <- sub(".*\\|", "", diaWorkflowResults[["Skyline_PROSIT_EDIA_GPF"]]$Protein.Names)
diaWorkflowResults[["Skyline_DIANN_AI_GPF"]]$Protein.Names <- sub(".*\\|", "", diaWorkflowResults[["Skyline_DIANN_AI_GPF"]]$Protein.Names)
diaWorkflowResults[["Spectronaut_MSFragger"]]$Protein.Names <- gsub("'", '', sub(".*\\|", "", diaWorkflowResults[["Spectronaut_MSFragger"]]$Protein.Names))
diaWorkflowResults[["Skyline_MSFragger"]]$Protein.Names <- sub(".*\\|", "", diaWorkflowResults[["Skyline_MSFragger"]]$Protein.Names)

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

  df <- df[,grepl("1-6", colnames(df)) | grepl("1-12", colnames(df)) | grepl("1-25", colnames(df))]
  return(df)
}

diaWorkflowResults <- lapply(diaWorkflowResults, function(x) {

  # remove proteins with both human and e.coli annotation
  x <- x[!(grepl("HUMAN", x$Protein.Names) & grepl("ECOLI", x$Protein.Names)),]

  x <- x[with(x, order(Protein.Names, Sample)),]
  x$Species <- NULL
  # If there are duplicates in terms of Protein.Names and Samples combined, check if one them is NA for PG.Quantity and remove respective entry
  # 9, 10, 11, 16, 17
  if (sum(duplicated(x[,c('Protein.Names', 'Sample')]))>0){
    frequencies <- x %>%
      group_by(Protein.Names, Sample) %>%
      summarise(Freq = n())

    x.freqAdded <- merge(x, frequencies)
    noDuplicate <- x.freqAdded$Freq == 1
    missing <- is.na(x$PG.Quantity)
    x <- x[noDuplicate | (!noDuplicate & !missing), ]
  }
  x <- generateWideFormattedDf(x)

  x
})

cols <- colnames(diaWorkflowResults[[1]])
design <- data.frame(colnames=cols, dilution=cols)

#design[!(grepl("1-6", design$colnames) | grepl("1-12", design$colnames) | grepl("1-25", design$colnames)), ]$dilution <- "undiluted"
design[grepl("1-6", design$colnames), ]$dilution <- "X6"
design[grepl("1-12", design$colnames), ]$dilution <- "X12"
design[grepl("1-25", design$colnames), ]$dilution <- "X25"
row.names(design) <- design$colnames

design$dilution <- factor(design$dilution,
                          levels = c("X25", "X12", "X6"))
design <- design[order(design$dilution), ]
design$colnames <- NULL

sort.idx <- row.names(design)
idx25 <- which(design$dilution == "X25")
idx12 <- which(design$dilution == "X12")

diaWorkflowResults <- lapply(diaWorkflowResults,function(x, sort.idx, idx25, idx12) {
  # Keep only samples of dilution 12 or 25
  x <- x[, c(idx25, idx12)]

  # Remove empty rows
  x <- x[rowSums(is.na(x)) != ncol(x),]

  sort.idx <- sort.idx[!grepl("1-6", sort.idx)]
  x <- x[,sort.idx]
  x
}, sort.idx=sort.idx, idx25=idx25, idx12 = idx12)

set.seed(50)
indices <- list()
# 100 iterations
for (i in seq(1:100)){
  # Sample sizes from 3 to 23
  for(sampleSize in seq(3, length(idx25))){
    indices25 <- sample(idx25, sampleSize, replace = T)
    indices12 <- sample(idx12, sampleSize, replace = T)
    indices <- list.append(indices, c(indices25, indices12))
  }
}

# same sample indices are used for all DIA workflows
repList <- list()
repList <- lapply(diaWorkflowResults,function(x, repList, indices) {
  lst <- lapply(indices,function(idx, x) {
    x <- x[,idx]
    x
  },  x=x)

  repList <- list.append(repList, lst)
  repList
}, repList=repList, indices=indices)


repList <- unlist(repList, recursive = FALSE)

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
save(diaWorkflowResults, indices, intersectProteinNames, combinedProteinNames, file="DIAsoftwareOutputProteinLevel_1to12And1to25Only_wideFormat_withBootstrapIndicesAndIntersectAndCombinedProteinNames.RData")

saveRDS(intersectProteinNames, file = "intersectProteinNames.rds")
saveRDS(combinedProteinNames, file = "combinedProteinNames.rds")
saveRDS(indices, file = "indices.rds")

session <- sessionInfo()
sink(paste0("sessionInfo_repList.txt"))
print(session)
sink()

# Produces 17 files with approx. 2-3 GB of size each
for (dianame in names(repList)){
  saveRDS(repList[[dianame]], file = paste0(dianame, ".rds"))
}
