# library(data.table)
library(rlist)
library(openxlsx)
library(svglite)
library(readr) # Dev version: devtools::install_github("tidyverse/readr")
library(forcats)

library(tidyr)
library(plyr)
library(dplyr)
library(reshape2)

library(ggplot2)
theme_set(theme_bw())
library(RColorBrewer)
library(cowplot)
library(gridExtra)
library(patchwork)
library(viridis)
library(ggridges)
library(UpSetR)
# library(grid)
# library(GGally)
library(corrplot)
library(corrr)

library(pcaMethods)
library(limma)

library(foreach)
library(doParallel)

my_greens = brewer.pal(n = 9, "Greens")[3:9]

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

addDilutionColumn <- function(df, colname) {
  df$dilution <- ""
  df[!(grepl("1-6", df[[colname]]) | grepl("1-12", df[[colname]]) | grepl("1-25", df[[colname]])), ]$dilution <- "Lymph nodes"
  df[grepl("1-6", df[[colname]]), ]$dilution <- "Lymph nodes + 1:6 E.coli"
  df[grepl("1-12", df[[colname]]), ]$dilution <- "Lymph nodes + 1:12 E.coli"
  df[grepl("1-25", df[[colname]]), ]$dilution <- "Lymph nodes + 1:25 E.coli"
  df$dilution <- factor(df$dilution, 
                        levels = c("Lymph nodes", "Lymph nodes + 1:25 E.coli", "Lymph nodes + 1:12 E.coli", "Lymph nodes + 1:6 E.coli" ))
  df
}

addSoftwareLibraryColumns <- function(df.combined.LibSoftware, sepSpace=FALSE) {
  if (sepSpace){
    df.combined.LibSoftware$dia <- gsub(" ", "_", df.combined.LibSoftware$dia)
  }
  
  uniqueDias <- unique(df.combined.LibSoftware$dia)
  # Replace DIANN_DIANN_AI by DIANN_Predicted and Spectronaut_DirectDIA by Spectronaut_Predicted
  if ("DIANN_DIANN_AI" %in% uniqueDias){
    df.combined.LibSoftware[df.combined.LibSoftware$dia == "DIANN_DIANN_AI",]$dia <- "DIANN_Predicted"
  } 
  if ("Spectronaut_DirectDIA" %in% uniqueDias){
    df.combined.LibSoftware[df.combined.LibSoftware$dia == "Spectronaut_DirectDIA",]$dia <- "Spectronaut_Predicted"
  }
  
  df.combined.LibSoftware$diaSoftware <- sub("_.*", "", df.combined.LibSoftware$dia)
  df.combined.LibSoftware$diaLibrary <- sub("^[^_]*_", "", df.combined.LibSoftware$dia)
  
  df.combined.LibSoftware$diaSoftware  <- factor(df.combined.LibSoftware$diaSoftware , 
                                                 levels=c("DIANN", "Skyline", "OSW", "Spectronaut"))
  
  df.combined.LibSoftware$diaLibrary <- gsub("_", " ", df.combined.LibSoftware$diaLibrary)
  
  df.combined.LibSoftware$diaLibrary <- factor(df.combined.LibSoftware$diaLibrary, 
                                               levels = c("PROSIT EDIA GPF", "DIANN AI GPF", "MaxQuant", "MSFragger", "Predicted" ))
  if (sepSpace){
    df.combined.LibSoftware$dia <- gsub("_", " ", df.combined.LibSoftware$dia)
  }
  df.combined.LibSoftware
}

########################################################################

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
  
  design <- design[order(design$dilution, design$colnames), ]
  design$colnames <- NULL
  design
}

# Save to xlsx file, each tab corresponding to one DIA software-library combination
addToExcelWb <- function(out_xlsx, df_for_workbook, worksheetName, rowNames=FALSE){
  # Check to see if file doesn't exist
  if (!file.exists(out_xlsx))  {
    # Create workbook using openxlsx
    wb <- openxlsx::createWorkbook()
  } else {
    wb <- openxlsx::loadWorkbook(out_xlsx)
  }
  # Add worksheet
  openxlsx::addWorksheet(wb, worksheetName)
  # Write data frame to new worksheet
  openxlsx::writeData(wb, worksheetName, df_for_workbook, rowNames = rowNames, colNames = TRUE, keepNA = TRUE, na.string="NA")
  # Save file
  openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
}

#for (diaworkflow.name in names(diaWorkflowResults)){
#  addToExcelWb(paste0("diaWorkflowResults_log.xlsx"), diaWorkflowResults[[diaworkflow.name]], strtrim(diaworkflow.name, 31), rowNames = TRUE)
#}

#saveRDS(diaWorkflowResults, file = "diaWorkflowResults_log.rds")
########################################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

diaWorkflowResults.precursors <- unlist(lapply("DIAsoftwareOutputPrecursorLevel.RData", function(x) mget(load(x))), recursive = FALSE)
diaWorkflowResults.precursors <- rlist::list.remove(diaWorkflowResults.precursors, c(".Random.seed", "OSW_DIANN", "anotation_df2"))

oldNames <- c("DIANN_PROSIT", "DIANN_predicted", "Spectronaut_PROSIT", "DIANN_Maxquant")
newNames <- c("DIANN_PROSIT_EDIA_GPF", "DIANN_DIANN_AI", "Spectronaut_PROSIT_EDIA_GPF", "DIANN_MaxQuant")

mm <- match(names(diaWorkflowResults.precursors), oldNames)
names(diaWorkflowResults.precursors)[!is.na(mm)] <- as.character(newNames[na.omit(mm)])

diaWorkflowResults.precursors <- diaWorkflowResults.precursors[order(names(diaWorkflowResults.precursors))]

formatMSTable.precursors <- function(x) {
  x <- as.data.frame(x)
  x <- x[, c("Precursor.ID", "Run", "Precursor.Quant", "Species")]
  colnames(x) <- c("Precursor.ID", "Sample", "Precursor.Quant", "Species")
  
  if (startsWith(toString(x[1,]$Sample), "LE_")){
    x$Sample <- gsub("LE_", "", x$Sample)
    x$Sample <- sub("_[^_]+$", "", x$Sample)
  }
  
  x$Sample <- gsub("D:\\\\data\\\\Klemens\\\\RealDilutionSeries_mzml\\\\Lymph_|D:\\\\data\\\\Klemens\\\\RealDilutionSeries_dia\\\\Lymph_", "", x$Sample)

  x$Sample <- gsub("in\\/Lymph_Ecoli_|in\\/Lymph_", "", x$Sample)
  x$Sample <- gsub("Lymph_|Lymph_Ecoli_|Ecoli_|Lymphnode_|_mzML\\.|\\.mzML|mzML|\\.dia", "", x$Sample)
  x$Sample <- gsub("0100", "100", x$Sample)
  
  x <- x[!(x$Species == "synthetic" | x$Species == "usion" | is.na(x$Species)),]
  
  # Add species information to Precursor.ID so that for plotting human and E. coli precursors can be separated
  x$Precursor.ID <- paste0(x$Precursor.ID, "_", x$Species)
  x
}


#diaWorkflowResults.precursors <- lapply(diaWorkflowResults.precursors, formatMSTable.precursors)
diaWorkflowResults.precursors <- parallel::mclapply(diaWorkflowResults.precursors, formatMSTable.precursors, mc.cores = detectCores() - 1)


# addLeadingZeros <- function(el){
#   if (length(el) == 1){
#     elFinal <- formatC(as.numeric(el),width=3,format='f',digits=0,flag='0')
#   } else if (length(el) == 2){
#     elZerosAdded <- formatC(as.numeric(el[2]),width=3,format='f',digits=0,flag='0')
#     elFinal <- paste0(el[1], "_", elZerosAdded)
#   }
#   elFinal
# }

# for (dia in c("OSW_MaxQuant", "OSW_DIANN_AI_GPF", "OSW_MSFragger")){
#   diaWorkflowResults.precursors[[dia]]$Sample <- 
#     unlist(lapply( strsplit(
#       unique(diaWorkflowResults.precursors[[dia]]$Sample),split='_', fixed=TRUE), addLeadingZeros))
#   
# }

# diaWorkflowResults.precursors[["OSW_MaxQuant"]]$Sample <- unlist(lapply( strsplit(unique(diaWorkflowResults.precursors[["OSW_MaxQuant"]]$Sample),split='_', fixed=TRUE), addLeadingZeros))
# diaWorkflowResults.precursors[["OSW_DIANN_AI_GPF"]]$Sample <- unlist(lapply( strsplit(unique(diaWorkflowResults.precursors[["OSW_DIANN_AI_GPF"]]$Sample),split='_', fixed=TRUE), addLeadingZeros))
# diaWorkflowResults.precursors[["OSW_MSFragger"]]$Sample <- unlist(lapply( strsplit(unique(diaWorkflowResults.precursors[["OSW_MSFragger"]]$Sample),split='_', fixed=TRUE), addLeadingZeros))


generateWideFormattedDf.precursors <- function(df) {
  # remove duplicates from multiple precursors carrying the same quantitative information for a certain protein
  df <- dplyr::distinct(df)
  
  df <- reshape(df, idvar="Precursor.ID", timevar="Sample",  direction="wide")
  colnames(df) <- sub("Precursor.Quant.", "", colnames(df))
  
  row.names(df) <- df$Precursor.ID
  df$Precursor.ID <- NULL
  df[df == 0] <- NA
  df <- log2(df)
  
  if (!"1-6_028" %in% names(df)){
    df["1-6_028"] <- NA
  }
  
  return(df)
}


cleanData.precursors <- function(x) {
  # remove proteins with both human and e.coli annotation
  x <- x[!(grepl("HUMAN", x$Precursor.ID) & grepl("ECOLI", x$Precursor.ID)),]
  
  x <- x[with(x, order(Precursor.ID, Sample)),]
  x$Species <- NULL
  # If there are duplicates in terms of Precursor.ID and Samples combined, check if one them is NA for PG.Quantity and remove respective entry
  if (sum(duplicated(x[,c('Precursor.ID', 'Sample')]))>0){
    frequencies <- x %>%
      group_by(Precursor.ID, Sample) %>%
      dplyr::summarise(Freq = n())
    
    x.freqAdded <- merge(x, frequencies)
    noDuplicate <- x.freqAdded$Freq == 1
    missing <- is.na(x$Precursor.Quant)
    x <- x[noDuplicate | (!noDuplicate & !missing), ]
  }
  x <- generateWideFormattedDf.precursors(x)
  
  # Remove empty rows
  x <- x[rowSums(is.na(x)) != ncol(x),]
  x
}

#diaWorkflowResults.precursors <- lapply(diaWorkflowResults.precursors, cleanData.precursors)

diaWorkflowResults.precursors <- parallel::mclapply(diaWorkflowResults.precursors, cleanData.precursors, mc.cores = detectCores() - 1)


design <- getDesign(diaWorkflowResults.precursors[[1]])

samplenames <- row.names(design)

diaWorkflowResults.precursors <- lapply(diaWorkflowResults.precursors,function(x, samplenames, idx25, idx12) {
  # Bring to same order as design data frame
  x <- x[, samplenames]
  x
}, samplenames=samplenames)

saveRDS(diaWorkflowResults.precursors, file = "diaWorkflowResults_allDilutions_precursorLevel.rds")

########################################################################
getStackedList <- function(idx, colID, diaWorkflowResults){
  x <- diaWorkflowResults[[idx]]
  x[colID] <- row.names(x)
  x$Species <- "Human"
  x[grepl("_ECOLI", x[[colID]]),]$Species <- "E. coli"
  x <- reshape2::melt(x)
  x <- addDilutionColumn(x, "variable") 
  
  #x$Species <- "Human"
  #x[grepl("_ECOLI", x[colID]),]$Species <- "E. coli"
  
  x$dia <- names(diaWorkflowResults)[idx]
  colnames(x) <- c(colID, "species", "sample", "intensity", "dilution", "dia")
  x$dia <- gsub("_", " ", x$dia )
  x <- addSoftwareLibraryColumns(x, sepSpace = TRUE)
  x$dilution <- factor(x$dilution,
                       levels = c("Lymph nodes", "Lymph nodes + 1:25 E.coli", "Lymph nodes + 1:12 E.coli", "Lymph nodes + 1:6 E.coli"))
  
  x
}


plotDIAWorkflowDistributions <- function(diaWorkflowResults, level=c("proteinLevel", "precursorLevel"), fileprefix="") {
  if (level == "proteinLevel"){
    colID <- "ProteinName"
  } else if (level == "precursorLevel"){
    colID <- "Precursor.ID"
  }

  stacked1 <- getStackedList(1, colID, diaWorkflowResults)
  stackedRest <- parallel::mclapply(2:length(diaWorkflowResults), getStackedList, colID = colID, diaWorkflowResults = diaWorkflowResults, mc.cores=detectCores()-1)
  stacked.df <- dplyr::bind_rows(append(stackedRest, list(stacked1)))

  # stacked.df$dia <- gsub("_", " ", stacked.df$dia )
  # stacked.df <- addSoftwareLibraryColumns(stacked.df, sepSpace = TRUE)
  # stacked.df$dilution <- factor(stacked.df$dilution,
  #                               levels = c("Lymph nodes", "Lymph nodes + 1:25 E.coli", "Lymph nodes + 1:12 E.coli", "Lymph nodes + 1:6 E.coli"))
  # 
  #addToExcelWb("diaWorkflowResults_log_longFormat.xlsx", stacked.df, strtrim("diaWorkflowResults_log_longFormat", 31), rowNames = FALSE)
  
  stacked.nonNA <- stacked.df %>% filter(!is.na(intensity)) %>% dplyr::group_by(sample, dilution, dia) %>% 
    dplyr::summarise(NumberProteins = n())
  stacked.nonNA.meanMedians <- stacked.nonNA %>%  dplyr::group_by(dilution, dia) %>% 
    dplyr::summarise(median = median(NumberProteins), mean = mean(NumberProteins))
  
  stacked.nonNA.speciesSeparated <- stacked.df %>% filter(!is.na(intensity)) %>% dplyr::group_by(sample, dilution, dia, species) %>% 
    dplyr::summarise(NumberProteins = n())
  
  if (level == "proteinLevel"){
    axisLabel <- "# Proteins"
  } else if (level =="precursorLevel"){
    axisLabel <- "# Precursors"
  }
  
  # FIGURE 2
  ggplot.intensity.violinBoxplot <- ggplot(stacked.nonNA,aes(forcats::fct_rev(dia), NumberProteins, fill=forcats::fct_rev(dilution))) +
    geom_boxplot(alpha=0.5, outlier.size=0.5) +
    coord_flip() +
    xlab("") +
    ylab(axisLabel) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_brewer(palette="Greens", direction=-1) +
    theme(legend.title = element_blank(), axis.text.y = element_text(hjust=0), legend.position="bottom",
          legend.justification = "left", legend.direction = "vertical") +
    guides(fill = guide_legend(reverse = TRUE)) 
  ggsave(file=paste0(fileprefix, "_ProteinNumbers_", level, ".svg"), ggplot.intensity.violinBoxplot, height=8, width=6)
  
  ggplot.intensity.violinBoxplot.Ecoli <- ggplot(stacked.nonNA.speciesSeparated[stacked.nonNA.speciesSeparated$species=="E. coli",],aes(forcats::fct_rev(dia), NumberProteins, fill=forcats::fct_rev(dilution))) +
    geom_boxplot(alpha=0.5, outlier.size=0.5) +
    coord_flip() +
    xlab("") +
    ylab(axisLabel) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_brewer(palette="Greens", direction=-1) +
    theme(legend.title = element_blank(), axis.text.y = element_text(hjust=0), legend.position="bottom",
          legend.justification = "left", legend.direction = "vertical") +
    guides(fill = guide_legend(reverse = TRUE)) 
  ggsave(file=paste0(fileprefix, "_ProteinNumbers_Ecoli_", level, ".svg"), ggplot.intensity.violinBoxplot.Ecoli, height=8, width=6)
  
  ggplot.intensity.violinBoxplot.Human <- ggplot(stacked.nonNA.speciesSeparated[stacked.nonNA.speciesSeparated$species=="Human",],aes(forcats::fct_rev(dia), NumberProteins, fill=forcats::fct_rev(dilution))) +
    geom_boxplot(alpha=0.5, outlier.size=0.5) +
    coord_flip() +
    xlab("") +
    ylab(axisLabel) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_brewer(palette="Greens", direction=-1) +
    theme(legend.title = element_blank(), axis.text.y = element_text(hjust=0), legend.position="bottom",
          legend.justification = "left", legend.direction = "vertical") +
    guides(fill = guide_legend(reverse = TRUE)) 
  ggsave(file=paste0(fileprefix, "_ProteinNumbers_Human_", level, ".svg"), ggplot.intensity.violinBoxplot.Human, height=8, width=6)
  
  
  ggplot.intensityFacet <- ggplot(stacked.df, aes(x = intensity, y=forcats::fct_rev(dia), fill=forcats::fct_rev(dilution))) + 
    ggridges::geom_density_ridges(alpha=0.5) +
    facet_grid(. ~ species) +
    ggplot2::xlab("Log2(Intensity)") +
    ggplot2::ylab("") +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_brewer(palette="Greens", direction=-1) +
    theme(legend.title = element_blank(), axis.text.y = element_text(hjust=0), legend.position="bottom",
          legend.justification = "left", legend.direction = "vertical") +
    guides(fill = guide_legend(reverse = TRUE)) 
  ggsave(file=paste0(fileprefix, "_intensityDistributions_", level, ".svg"), ggplot.intensityFacet, height=8, width=6)
  
  # pdf(file="diaWorkflowDensityPlot_grid_Ecoli.pdf", width=12, height=8)
  # ggplot(stacked.df[stacked.df$species == "E. coli", ], aes(x =  intensity, fill=dilution)) + 
  #   geom_density(aes(fill=dilution), alpha = 0.5) + 
  #   # facet_wrap(~dia) +
  #   facet_grid(diaSoftware ~ diaLibrary, scales = "free") + 
  #   ggplot2::scale_fill_brewer(palette=rev("Greens")) +
  #   theme(legend.title = element_blank()) + 
  #   xlab("Log2(Intensity)") +
  #   ylab("Density") +
  #   ggtitle("E. coli")
  # dev.off()
  # 
  # pdf(file="diaWorkflowDensityPlot_grid_Human.pdf", width=12, height=8)
  # ggplot(stacked.df[stacked.df$species == "Human", ], aes(x =  intensity, fill=dilution)) + 
  #   geom_density(aes(fill=dilution), alpha = 0.5) + 
  #   # facet_wrap(~dia) +
  #   facet_grid(diaSoftware ~ diaLibrary, scales = "free") + 
  #   ggplot2::scale_fill_brewer(palette=rev("Greens")) +
  #   theme(legend.title = element_blank()) + 
  #   xlab("Log2(Intensity)") +
  #   ylab("Density") +
  #   ggtitle("Human")
  # dev.off()
  
  stacked.df.var <- stacked.df %>% 
    group_by(!!as.name(colID), dilution, dia, diaSoftware, diaLibrary, species) %>% 
    dplyr::summarize(variance = var(intensity, na.rm = TRUE))
  
  ggplot.var <- ggplot(stacked.df.var[stacked.df.var$species == "E. coli", ], aes(x = log2(variance), 
                                                                                  y=forcats::fct_rev(dia), fill=forcats::fct_rev(dilution))) + 
    geom_boxplot(alpha=0.5, outlier.size=0.5) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_brewer(palette="Greens", direction=-1) +
    theme(legend.title = element_blank(), axis.text.y = element_text(hjust=0), legend.position="bottom",
          legend.justification = "left", legend.direction = "vertical") +
    guides(fill = guide_legend(reverse = TRUE)) +
    xlab("Log2(Variance)") +
    ylab("") +
    ggtitle("E. coli")
  ggsave(file=paste0("proteinVariance_Ecoli_", level, ".svg"), ggplot.var, height=8, width=6)
  
  # FIGURE 2
  ggplot.var.mins15 <- ggplot.var + xlim(-15, NA) 
  ggsave(file=paste0(fileprefix, "_proteinVariance_Ecoli_xlimMinus15_", level, ".svg"), ggplot.var.mins15, height=8, width=6)
  
  
  ggplot.var.human <- ggplot(stacked.df.var[stacked.df.var$species == "Human", ], aes(x = log2(variance), 
                                                                                      y=forcats::fct_rev(dia), fill=forcats::fct_rev(dilution))) + 
    geom_boxplot(alpha=0.5, outlier.size=0.5) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_brewer(palette="Greens", direction=-1) +
    theme(legend.title = element_blank(), axis.text.y = element_text(hjust=0), legend.position="bottom",
          legend.justification = "left", legend.direction = "vertical") +
    guides(fill = guide_legend(reverse = TRUE)) +
    xlab("Log2(Variance)") +
    ylab("") +
    ggtitle("Human")
  ggsave(file=paste0("proteinVariance_Human_", level, ".svg"), ggplot.var.human, height=8, width=6)
  
  ggplot.var.human.mins20 <- ggplot.var.human + xlim(-20, NA) 
  ggsave(file=paste0(fileprefix, "_proteinVariance_Human_xlimMinus20_", level, ".svg"), ggplot.var.human.mins20, height=8, width=6)
}
########################################################################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

diaWorkflowResults <- readRDS(file = "diaWorkflowResults_allDilutions.rds")

diaWorkflowResults <- diaWorkflowResults[!names(diaWorkflowResults) %in% c("SkylineDIANN_DIANN_AI_GPF", "SkylineDIANN_MaxQuant", 
                                                                           "SkylineDIANN_MSFragger", "SkylineDIANN_PROSIT_EDIA_GPF")]
names(diaWorkflowResults) <- sub("SkylineMSstats", "Skyline",names(diaWorkflowResults))
names(diaWorkflowResults)[names(diaWorkflowResults)=="DIANN_DIANN_AI"] <- "DIANN_Predicted"
names(diaWorkflowResults)[names(diaWorkflowResults)=="Spectronaut_DirectDIA"] <- "Spectronaut_Predicted"
diaWorkflowResults <- diaWorkflowResults[order(names(diaWorkflowResults))]


# Remove data of sample 1-6_028 due to too many missing values
diaWorkflowResults <- lapply(diaWorkflowResults, function(x){
  x$`1-6_028` <- NULL
  x
})

diaWorkflowResults.precursors <- lapply(diaWorkflowResults.precursors, function(x){
  x$`1-6_028` <- NULL
  x
})

plotDIAWorkflowDistributions(diaWorkflowResults, level = "proteinLevel", fileprefix="Fig2")
plotDIAWorkflowDistributions(diaWorkflowResults.precursors, level= "precursorLevel", fileprefix="FigS3")


########################################################################
# Process TRIC 0.05 and TRIC 0.01 of OpenSwath
tric001 <- read.csv("Proteinexpressionmatrix_MaxQuant_TRIC_0_01_max-fdr-quality.txt", check.names = FALSE, sep="\t")
tric005 <- read.csv("Proteinexpressionmatrix_MaxQuant_TRIC_0_05_max-fdr-quality.txt", check.names = FALSE, sep="\t")

addLeadingZeros <- function(el){
  if (length(el) == 1){
    elFinal <- formatC(as.numeric(el),width=3,format='f',digits=0,flag='0')
  } else if (length(el) == 2){
    elZerosAdded <- formatC(as.numeric(el[2]),width=3,format='f',digits=0,flag='0')
    elFinal <- paste0(el[1], "_", elZerosAdded)
  }
  elFinal
}

cleanTricDataframes <- function(expressionDf){
  row.names(expressionDf) <- expressionDf$ProteinName
  expressionDf$ProteinName <- NULL
  
  df_colnames <- colnames(expressionDf)
  
  if (startsWith(toString(df_colnames), "LE_")){
    df_colnames <- gsub("LE_", "", df_colnames)
    df_colnames <- sub("_[^_]+$", "", df_colnames)
  }
  df_colnames <- gsub("Lymph_|Lymph_Ecoli_|Ecoli_|Lymphnode_|\\.mzML|mzML", "", df_colnames)
  #x$Sample <- gsub("0100", "100", x$Sample)
  

  df_colnames <- unlist(lapply( strsplit(df_colnames,split='_', fixed=TRUE), addLeadingZeros))
  colnames(expressionDf) <- df_colnames
  expressionDf[expressionDf == 0] <- NA
  expressionDf <- log2(expressionDf)
  expressionDf$`1-6_028` <- NULL
  expressionDf <- expressionDf[rowSums(is.na(expressionDf)) != ncol(expressionDf),]
  expressionDf
}

tric001 <- cleanTricDataframes(tric001)
tric005 <- cleanTricDataframes(tric005)

tric_rownames_unique <- unique(c(row.names(tric001), row.names(tric005)))

# data.table::fwrite(list(tric_rownames_unique), file = "OSW_MaxQuant_tric_unique.txt")
## Derived translation table from https://www.uniprot.org/uploadlists/

mappedProteins <- read.table("OSW_MaxQuant_tric_unique_mappedProteins_20220203.tsv", sep = "\t", header = TRUE)
row.names(tric001) <- mappedProteins$Entry.name[match(row.names(tric001), mappedProteins$Entry)]
row.names(tric005) <- mappedProteins$Entry.name[match(row.names(tric005), mappedProteins$Entry)]


design <- getDesign(diaWorkflowResults[[1]])

samplenames <- row.names(design)
tric001 <- tric001[, samplenames]
tric005 <- tric005[, samplenames]

########################################################################
getProteinMeansNAsForDf <- function(df, dia="") {
  proteinMeans <- rowMeans(df, na.rm=TRUE)
  percNAsInProtein <- (rowSums(is.na(df))/ncol(df)) * 100
  protein.df <- data.frame(percNAsInProtein, proteinMeans, dia, proteinName=row.names(df))
  protein.df
}

getProteinMeansNAs <- function(idx, diaWorkflowResults){
  dia <- names(diaWorkflowResults)[idx]
  print(dia)
  df <- diaWorkflowResults[[idx]]
  protein.df <- getProteinMeansNAsForDf(df, dia) 
  protein.df
}

getSampleNADf <-  function(idx, diaWorkflowResults, species=""){
  dia <- names(diaWorkflowResults)[idx]
  print(dia)
  df <- diaWorkflowResults[[idx]]
  
  if (species != ""){
    df <- df[grepl(paste0("_", species), row.names(df)),]
  }
  
  sampleMeans <- colMeans(df, na.rm=TRUE)
  percNAsInSample <- (colSums(is.na(df))/nrow(df)) * 100
  sample.df <- data.frame(percNAsInSample, sampleMeans, dia, sampleName=colnames(df))
  
  sample.df
}

getSampleNADfPerSpecies <- function(diaWorkflowResults, species="") {
  sampleMeansNAsPlot.humanAndEcoli.df <- do.call(rbind, 
                                                 lapply(1:length(diaWorkflowResults), 
                                                        getSampleNADf, 
                                                        diaWorkflowResults=diaWorkflowResults, 
                                                        species=species))
  sampleMeansNAsPlot.humanAndEcoli.df <- addDilutionColumn(sampleMeansNAsPlot.humanAndEcoli.df, "sampleName") 
  sampleMeansNAsPlot.humanAndEcoli.df
}

plotMeanCorrelations <- function(diaWorkflowResults, level="", prefix="") {
  proteinMeansNAsPlot.lst <- lapply(1:length(diaWorkflowResults), getProteinMeansNAs, diaWorkflowResults=diaWorkflowResults)
  
  proteinMeansNA.df <- do.call(rbind, proteinMeansNAsPlot.lst)
  proteinMeansNA.df$species <- "Human"
  proteinMeansNA.df[grepl("_ECOLI", proteinMeansNA.df$proteinName),]$species <- "E. coli"
  
  for (species in c("E. coli", "Human", "E. coli & Human")) {
    df <- proteinMeansNA.df
    if (species != "E. coli & Human"){
      df <- df[df$species==species,]
    }
    
    df <- addSoftwareLibraryColumns(df, sepSpace = TRUE) 
    
    # FIGURE 3A
    pdf(file=paste0(prefix, "A_proteinMeansNACorrelation_", gsub(" |\\.", "_", species), "_", level, ".pdf"), width=8, height=6)
    if (level=="precursorLevel"){
      ylab.name <- "Precursor mean"
    } else {
      ylab.name <- "Protein mean"
    }
    
    if (species == "E. coli & Human"){
      if (level=="precursorLevel"){
        alpha <- .002
      } else {
        alpha <- .05 
      }
      print(ggplot(df, aes(x = percNAsInProtein, y = proteinMeans, color=forcats::fct_rev(species))) + 
              geom_point(size=0.5, alpha=alpha) +
              facet_grid(diaSoftware ~ diaLibrary, scales = "free") +
              scale_color_manual(values=c(brewer.pal(11, "RdYlBu")[10], brewer.pal(11, "RdYlBu")[2])) +
              xlab("% Missing values") +
              ylab(ylab.name) +
              ggtitle(paste0(species, ", ", level)) +
              theme(legend.title = element_blank(), axis.text.y = element_text(hjust=0), legend.position="bottom",
                    legend.justification = "left", legend.direction = "vertical") +
              guides(colour = guide_legend(override.aes = list(alpha = 1))))
      
    } else {
      print(ggplot(df, aes(x = percNAsInProtein, y = proteinMeans)) + 
              geom_point(size=0.5, alpha=.05) +
              facet_grid(diaSoftware ~ diaLibrary, scales = "free") +
              xlab("% Missing values") +
              ylab(ylab.name) +
              ggtitle(paste0(species, ", ", level)) +
              theme(legend.title = element_blank(), axis.text.y = element_text(hjust=0), legend.position="bottom",
                    legend.justification = "left", legend.direction = "vertical") )
    }
    dev.off()
    
    if (species == "E. coli & Human"){
      speciesString <- ""
    } else if (species == "E. coli"){
      speciesString <- "ECOLI"
    } else if (species == "Human"){
      speciesString <- "HUMAN"
    }
    
    df2 <- getSampleNADfPerSpecies(diaWorkflowResults, species=speciesString)
    df2 <- addSoftwareLibraryColumns(df2, sepSpace = TRUE) 
    
    pdf(file=paste0(prefix, "B_sampleMeansNACorrelation_", gsub(" |\\.", "_", species), "_", level, ".pdf"), width=8, height=6)
    # FIGURE 3B
    print(ggplot(df2, aes(x = percNAsInSample, y = sampleMeans, color=forcats::fct_rev(dilution))) + 
            geom_point(size=0.9, alpha=0.5, stroke = 0.1) +
            facet_grid(diaSoftware ~ diaLibrary,  scales = "free") +
            xlab("% Missing values") +
            ylab("Sample mean") +
            ggtitle(paste0(species, ", ", level)) +
            theme(legend.title = element_blank(), axis.text.y = element_text(hjust=0), legend.position="bottom",
                  legend.justification = "left", legend.direction = "vertical") +
            guides(color = guide_legend(override.aes = list(alpha = 1))) )
    dev.off()
  }
}


plotMeanCorrelations(diaWorkflowResults, level="proteinLevel", prefix="Fig3") 
plotMeanCorrelations(diaWorkflowResults.precursors, level="precursorLevel", prefix="FigS7_FigS8_") 

##### Plotting TRIC
proteinMeansNAs.tric001 <- getProteinMeansNAsForDf(tric001, dia= "0.01 TRIC")
proteinMeansNAs.tric005 <- getProteinMeansNAsForDf(tric005, dia= "0.05 TRIC")
proteinMeansNAs.tricNo <- getProteinMeansNAsForDf(diaWorkflowResults[["OSW_MaxQuant"]], dia= "No TRIC")
proteinMeansNAs.tric.combined <- do.call("rbind", list(proteinMeansNAs.tric001, 
                                                      proteinMeansNAs.tric005,
                                                      proteinMeansNAs.tricNo))
proteinMeansNAs.tric.combined$species <- "Human"
proteinMeansNAs.tric.combined[grepl("_ECOLI", proteinMeansNAs.tric.combined$proteinName),]$species <- "E. coli"

svg(file=paste0("FigS6_TRIC_comparison_OSWMaxQuant.svg"), width=10, height=5)
# FIGURE 3B
print(ggplot(data = proteinMeansNAs.tric.combined, aes(x = percNAsInProtein, y = proteinMeans, color=forcats::fct_rev(species))) +
        geom_point(alpha=.1) +
        scale_color_manual(values=c(brewer.pal(11, "RdYlBu")[10], brewer.pal(11, "RdYlBu")[2])) +
        facet_wrap(~ dia) +
        xlab("% Missing values") +
        ylab("Protein mean") +
        theme(legend.title = element_blank(), axis.text.y = element_text(hjust=0), legend.position="bottom",
              legend.justification = "left", legend.direction = "vertical") +
        guides(colour = guide_legend(override.aes = list(alpha = 1)))
)
dev.off()
########################################################################

diaWorkflowResults.rbind <- do.call(rbind, diaWorkflowResults)

mtxTwoGroups <- diaWorkflowResults.rbind[,grepl("1-25|1-12", colnames(diaWorkflowResults.rbind))]
mtxTwoGroups <- mtxTwoGroups[rowSums(is.na(mtxTwoGroups)) != ncol(mtxTwoGroups), ]
group.size <- ncol(mtxTwoGroups)/2

log2FC.mtxTwoGroups <- data.frame(Protein=sub(".*\\.", "", row.names(mtxTwoGroups)), dia = sub("\\..*", "",  row.names(mtxTwoGroups)) , log2FC=rowMeans(mtxTwoGroups[, (group.size+1):ncol(mtxTwoGroups)], na.rm = TRUE)-rowMeans(mtxTwoGroups[, 1:group.size], na.rm = TRUE))
log2FC.mtxTwoGroups <- log2FC.mtxTwoGroups[!is.na(log2FC.mtxTwoGroups$log2FC),]
log2FC.mtxTwoGroups <- addSoftwareLibraryColumns(log2FC.mtxTwoGroups, sepSpace = TRUE) 

df.ecoli <- log2FC.mtxTwoGroups[grepl("ECOLI", log2FC.mtxTwoGroups$Protein), ]
df.human <- log2FC.mtxTwoGroups[grepl("HUMAN", log2FC.mtxTwoGroups$Protein), ]


# Label mode
densMode <- function(x){
  td <- density(x, na.rm=TRUE)
  maxDens <- which.max(td$y)
  list(x=td$x[maxDens], y=td$y[maxDens])
}

df.ecoli.mode <- ddply(df.ecoli,"dia", mutate, val_mean = signif(densMode(log2FC)$x,3), med.x = signif(densMode(log2FC)$x,3), med.y=signif(densMode(log2FC)$y,3)*1.1)

# SUPPLEMENTARY FIGURE S12
svg(file="FigS12_FC_ecoli.svg", width=12, height=8)
ggplot(data=df.ecoli.mode, aes(x=log2FC)) + 
  geom_density() + 
  geom_vline(aes(xintercept=val_mean), df.ecoli.mode, color="red",linetype="dashed",size=0.5) + 
  facet_grid(diaSoftware ~ diaLibrary, scales = "free_y") + 
  geom_text(data = df.ecoli.mode, aes(x=med.x,y=med.y,label=val_mean)) +
  xlab("Log2(Fold change)") +
  ylab("Density")
dev.off()


# df.human.mode <- ddply(df.human,"dia", mutate, val_mean = signif(densMode(log2FC)$x,3), med.x = signif(densMode(log2FC)$x,3), med.y=signif(densMode(log2FC)$y,3)*1.1)
# pdf(file="diaWorkflowViolinPlot_FC_human.pdf", width=12, height=8)
# ggplot(data=df.human.mode, aes(x=log2FC)) + 
#   geom_density() + 
#   geom_vline(aes(xintercept=val_mean), df.human.mode, color="red",linetype="dashed",size=0.5) + 
#   facet_grid(diaSoftware ~ diaLibrary, scales = "free_y", drop=TRUE) + 
#   geom_text(data = df.human.mode, aes(x=med.x,y=med.y,label=val_mean)) +
#   xlab("Log2(Fold change)") +
#   ylab("Density")
# dev.off()

########################################################################


getIntensitiesPerDilution <- function(dataset, humanOrEcoli=""){
  #dilutions <- c("undiluted", "X25", "X12", "X6")
  dilutions <- c("Lymph nodes", "Lymph nodes + 1:25 E.coli", "Lymph nodes + 1:12 E.coli", "Lymph nodes + 1:6 E.coli" )
  design <- getDesign(dataset)
  
  dataset.intensityCombined.lst <- list()
  dataset.intensityProteinmean.lst <- list()
  for (dilution in dilutions) {
    #print(dilution)
    dil.index <- which(design$dilution == dilution)
    dil.df <- dataset[, dil.index]
    
    if (humanOrEcoli == "ECOLI") {
      dil.df <- dil.df[grepl("ECOLI", row.names(dil.df)),]
    } else if (humanOrEcoli == "HUMAN") {
      dil.df <- dil.df[grepl("HUMAN", row.names(dil.df)),]
    }
    
    dil.df.proteinMeans <- rowMeans(dil.df, na.rm = TRUE)
    dil.df.proteinMeans[is.nan(dil.df.proteinMeans)] <- NA
    dataset.intensityProteinmean.lst <- rlist::list.append(dataset.intensityProteinmean.lst, dil.df.proteinMeans)
    
    dil.df <- dil.df[rowSums(is.na(dil.df)) != ncol(dil.df),]
    dataset.intensityCombined.lst <- rlist::list.append(dataset.intensityCombined.lst, unlist(dil.df))
    
  }
  names(dataset.intensityCombined.lst) <- names(dataset.intensityProteinmean.lst) <- dilutions
  list(intensityCombined = dataset.intensityCombined.lst, intensityProteinMean = dataset.intensityProteinmean.lst)
}

# intensityCombOrMean: "intensityCombined" or "intensityProteinMean"
getCombIntensitiesMeans <- function(idx, intensitiesPerDilution, intensityCombOrMean) {
  x <- intensitiesPerDilution[[idx]]
  dia <- names(intensitiesPerDilution)[idx]
  if (intensityCombOrMean == "intensityCombined"){
    df <- stack(x[[intensityCombOrMean]])
  } else if (intensityCombOrMean == "intensityProteinMean"){
    df <- data.frame(x[["intensityProteinMean"]], check.names = FALSE)
  }
  
  df$dia <- dia
  df <- addSoftwareLibraryColumns(df, sepSpace = FALSE)
  df$diaSoftware <- as.character(df$diaSoftware)  
  df$diaLibrary <- as.character(df$diaLibrary)
  
  if (intensityCombOrMean == "intensityProteinMean"){
    df$ProteinNames <- row.names(df)
    df <-reshape2::melt(df, id.vars = c("ProteinNames","dia", "diaSoftware", "diaLibrary"))
  }
  df
}

getCombinedIntensitiesPerDIAWorkflow <- function(diaWorkflowResults, humanOrEcoli=""){
  intensitiesPerDilution <- parallel::mclapply(diaWorkflowResults, getIntensitiesPerDilution, humanOrEcoli=humanOrEcoli, mc.cores = detectCores() - 1)
  intensitiesPerDIAWorkflow <- parallel::mclapply(1:length(intensitiesPerDilution), getCombIntensitiesMeans, 
                                      intensitiesPerDilution=intensitiesPerDilution,intensityCombOrMean= "intensityCombined", mc.cores = detectCores() - 1)
  intensitiesPerDIAWorkflow.df <- do.call("rbind", intensitiesPerDIAWorkflow)
  colnames(intensitiesPerDIAWorkflow.df) <- c("intensity", "dilution", "dia", "diaSoftware", "diaLibrary")
  
  
  intensitiesPerDIAWorkflow.df$diaSoftware  <- factor(intensitiesPerDIAWorkflow.df$diaSoftware , 
                                                 levels=c("DIANN", "Skyline", "OSW", "Spectronaut"))
  
  intensitiesPerDIAWorkflow.df$diaLibrary <- factor(intensitiesPerDIAWorkflow.df$diaLibrary, 
                                               levels = c("PROSIT EDIA GPF", "DIANN AI GPF", "MaxQuant", "MSFragger", "Predicted" ))
  
  #intensitiesPerDIAWorkflow.df <- addSoftwareLibraryColumns(intensitiesPerDIAWorkflow.df, sepSpace = TRUE)
  intensitiesPerDIAWorkflow.df
}

getProteinMeansIntensitiesPerDIAWorkflow <- function(diaWorkflowResults, humanOrEcoli=""){
  proteinMeansPerDilution <- lapply(diaWorkflowResults, getIntensitiesPerDilution, humanOrEcoli=humanOrEcoli)
  proteinMeansPerDIAWorkflow <- lapply(1:length(proteinMeansPerDilution), getCombIntensitiesMeans, 
                                       intensitiesPerDilution=proteinMeansPerDilution, intensityCombOrMean= "intensityProteinMean")
  proteinMeansPerDIAWorkflow.df <- do.call("rbind", proteinMeansPerDIAWorkflow)
  colnames(proteinMeansPerDIAWorkflow.df) <- c("proteinnames", "dia", "diaSoftware", "diaLibrary", "dilution", "intensity")
  #proteinMeansPerDIAWorkflow.df <- addSoftwareLibraryColumns(proteinMeansPerDIAWorkflow.df, sepSpace = TRUE)
  proteinMeansPerDIAWorkflow.df
}

give.AverageNumNonNAsPerSample <- function(x){
  return(c(y = max(x)*1.07, label = round(sum(!is.na(x))/23, digits = 2)))  # 23 is the number of samples per dilution
}

plotIntensityViolinplotsInGrid <- function(diaWorkflowResults, level=c("proteinLevel", "precursorLevel"), fileprefix="") {
  comInt.Ecoli <- getCombinedIntensitiesPerDIAWorkflow(diaWorkflowResults, humanOrEcoli="ECOLI")
  comInt.Human <- getCombinedIntensitiesPerDIAWorkflow(diaWorkflowResults, humanOrEcoli="HUMAN")
  #comInt.EcoliAndHuman <- getCombinedIntensitiesPerDIAWorkflow(diaWorkflowResults, humanOrEcoli="")
  
  #protMeanInt.Ecoli <- getProteinMeansIntensitiesPerDIAWorkflow(diaWorkflowResults, humanOrEcoli="ECOLI")
  #protMeanInt.Human <- getProteinMeansIntensitiesPerDIAWorkflow(diaWorkflowResults, humanOrEcoli="HUMAN")
  #protMeanInt.EcoliAndHuman <- getProteinMeansIntensitiesPerDIAWorkflow(diaWorkflowResults, humanOrEcoli="")
  
  # pdf(file="diaWorkflowLinePlot_protMeanInt.Ecoli.pdf", width=10, height=8)
  # ggplot(protMeanInt.Ecoli, aes(x = dilution, y = intensity)) + 
  #   geom_line(aes(x = dilution, y = intensity, group=proteinnames), alpha=0.2, size=0.05) +
  #   facet_grid(diaSoftware ~ diaLibrary, scales = "fixed") +
  #   xlab("") +
  #   ylab("Log2(Intensity)") + 
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #   ggtitle("E.coli")
  # dev.off()
  # 
  # pdf(file="diaWorkflowLinePlot_protMeanInt.Human.pdf", width=10, height=8)
  # ggplot(protMeanInt.Human, aes(x = dilution, y = intensity)) + 
  #   geom_line(aes(x = dilution, y = intensity, group=proteinnames), alpha=0.2, size=0.05) +
  #   facet_grid(diaSoftware ~ diaLibrary, scales = "fixed") + 
  #   xlab("") +
  #   ylab("Log2(Intensity)") + 
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #   ggtitle("Human")
  # dev.off()
  
  # SUPPLEMENTARY FIGURE S4
  pdf(file=paste0(fileprefix, "_intensityViolinplots_forEachdiaWorkflow_Ecoli_", level,".pdf"), width = 10, height = 8)
  print(ggplot(comInt.Ecoli, aes(x = dilution, y = intensity)) + 
    geom_violin(scale="count") +
    geom_boxplot(width=0.08, outlier.size=0.5, outlier.alpha=0.4) +
    stat_summary(fun.data = give.AverageNumNonNAsPerSample, geom = "text", fun = max, size=2) +
    facet_grid(diaSoftware ~ diaLibrary, scales = "fixed") + 
    xlab("") +
    ylab("Log2(Intensity)") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_hline(aes(yintercept = median(intensity, na.rm = TRUE)), colour = 'red') +
    ggtitle(paste0("E.coli, ", level)))
  dev.off()
  
  pdf(file=paste0("S5_intensityViolinplots_forEachdiaWorkflow_human_", level, ".pdf"), width=10, height=8)
  print(ggplot(comInt.Human, aes(x = dilution, y = intensity)) + 
    geom_violin(scale="count") +
    geom_boxplot(width=0.08, outlier.size=0.5, outlier.alpha=0.4) +
    stat_summary(fun.data = give.AverageNumNonNAsPerSample, geom = "text", fun = max, size=2) +
    facet_grid(diaSoftware ~ diaLibrary, scales = "fixed") + 
    xlab("") +
    ylab("Log2(Intensity)") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_hline(aes(yintercept = median(intensity, na.rm = TRUE)), colour = 'red') +
    ggtitle(paste0("Human, ", level)))
  dev.off()
}


plotIntensityViolinplotsInGrid(diaWorkflowResults, level="proteinLevel", fileprefix="S4")
#plotIntensityViolinplotsInGrid(diaWorkflowResults.precursors, level="precursorLevel", fileprefix="S")

############################################################

getDfModel <- function(modus, df) {
  if (modus == "unnormalized"){
    df.model <- df
  } else if (modus == "TRQN") {
    mtx <- as.matrix(df)
    df.trqn <- MBQN::mbqn(mtx, FUN = mean)
    row.names(df.trqn) <- row.names(df)
    df.model <- as.data.frame(df.trqn)
  } else if (modus == "QN"){
    mtx <- as.matrix(df)
    df.qn <- MBQN::mbqn(mtx, FUN = NULL)
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

pdf(paste0("boxplotSpeciesSeparated.pdf"))
for (i in 1:length(diaWorkflowResults)) {
  dia.name <- names(diaWorkflowResults)[i]
  df <- diaWorkflowResults[[i]]
  df <- df[,grepl("1-25|1-12", colnames(df))]
  
  normalizations <- c("unnormalized", "TRQN", "QN", "median")
  
  for (normalization in normalizations){
    df.norm <- getDfModel(normalization, df)
    df.norm$species <- "Human"
    df.norm[grepl("_ECOLI", row.names(df.norm)),]$species <- "E. coli"
    df.norm.long <- reshape2::melt(df.norm)
    print(ggplot(df.norm.long, aes(x=value, y=forcats::fct_rev(variable), fill=species)) + 
            geom_boxplot() +
            geom_vline(aes(xintercept = median(value, na.rm = TRUE)), colour = 'red') +
            xlab("Log2(Intensity)") +
            ylab("DIA Workflow") +
            theme(legend.title = element_blank())+
            ggtitle(paste0(dia.name, ", ", normalization)))
  }
}
dev.off()


combinedProteinNames <- c()
combinedProteinNames <- lapply(diaWorkflowResults, function(x, combinedProteinNames){
  names <- unique(unlist(strsplit(row.names(x), ";")))
  combinedProteinNames <- c(combinedProteinNames, names)
  combinedProteinNames
}, combinedProteinNames=combinedProteinNames)

intersectProteinNames <- data.frame(table(unlist(combinedProteinNames)))
# Get protein names present in ALL DIA workflows
intersectProteinNamesAll <- intersectProteinNames[intersectProteinNames$Freq == length(diaWorkflowResults),]$Var1

old <- c("CTGEF_HUMAN;CTGE9_HUMAN;CTGE8_HUMAN;CTGE6_HUMAN;CTGE4_HUMAN;MIA2_HUMAN",
         "IMA6_HUMAN;IMA7_HUMAN",
         "LIN7A_HUMAN;LIN7C_HUMAN",
         "NDUCR_HUMAN;NDUC2_HUMAN",
         "NUD4B_HUMAN;NUDT4_HUMAN",
         "RRAGB_HUMAN;RRAGA_HUMAN",
         "SBNO1_HUMAN;SBNO2_HUMAN", 
         "WAC2A_HUMAN;WAC2C_HUMAN")
new <- c("MIA2_HUMAN",
         "IMA7_HUMAN",
         "LIN7C_HUMAN",
         "NDUC2_HUMAN",
         "NUDT4_HUMAN",
         "RRAGA_HUMAN",
         "SBNO1_HUMAN",
         "WAC2A_HUMAN")



diaWorkflowResults.intersect <- lapply(diaWorkflowResults, function(x){
  x <- x[row.names(x) !=  "LAP2A_HUMAN;LAP2B_HUMAN", ]
  for (j in 1:length(old)){
    #print(df$protein.name == old[j])
    try(row.names(x)[row.names(x) == old[j]] <- new[j])
  }
  
  x <- x[row.names(x) %in% intersectProteinNamesAll, ]
  x <- x[order(row.names(x)),]
  
  x
})

diaWorkflowResults.intersect <- lapply(diaWorkflowResults.intersect, function(x){
  unlist(x)
})

diaWorkflowResults.intersect.df <- do.call(cbind, diaWorkflowResults.intersect)  
colnames(diaWorkflowResults.intersect.df) <- names(diaWorkflowResults.intersect)


M <- cor(diaWorkflowResults.intersect.df, use="pairwise.complete.obs") # "pairwise.complete.obs" as there are NAs in regVals for some dias

svg(file = paste0("FigS2_corrplot_log2Intensities_diaworkflows.svg"), width = 16, height = 16)
corrplot::corrplot.mixed(M,  upper = "ellipse", lower = "number",
                         tl.pos = "lt", tl.col = "black", mar=c(0,0,2,0))
dev.off()

svg(file = paste0("networkplot_log2Intensities_diaworkflows.svg"), width = 25, height = 25)
diaWorkflowResults.intersect.df %>% 
  corrr::correlate() %>% 
  corrr::network_plot()
dev.off()


########################################################################

plotUpsetPlot <- function(diaWorkflowResults, selectedDilutions=c("Lymph nodes", "Lymph nodes + 1:25 E.coli", 
                                                                      "Lymph nodes + 1:12 E.coli", "Lymph nodes + 1:6 E.coli"), str="") {
  diaWorkflowResults.names <- lapply(diaWorkflowResults, function(df){
    design.dilutions <-  addDilutionColumn(data.frame(dia=colnames(df)), colname="dia") 
    
    df <- df[, design.dilutions$dilution %in% selectedDilutions]
    # Remove empty rows
    df <- df[rowSums(is.na(df)) != ncol(df),]
    names <- unique(unlist(strsplit(row.names(df), ";")))
    names
  })
  
  # SUPPLEMENTARY FIGURE S10 AND S11
  
  lapply(c("E.coli & Human", "E.coli", "Human"), function(intersectionType, diaWorkflowResults.names) {
    print(intersectionType)
    if (intersectionType == "E.coli"){
      protein.names <- lapply(diaWorkflowResults.names, function(x){
        x[grepl("ECOLI", x)]
      })
    } else if (intersectionType == "Human"){
      protein.names <- lapply(diaWorkflowResults.names, function(x){
        x[grepl("HUMAN", x)]
      })
    } else if (intersectionType == "E.coli & Human"){
      protein.names <- diaWorkflowResults.names
    }
    
    svg(paste0("FigS10AndS11_upsetPlot_", gsub(" ", "", intersectionType), "_noHEKIncluded_", str, ".svg"), width=10, height=8)
    print(UpSetR::upset(UpSetR::fromList(protein.names),
                        order.by = "freq", 
                        nsets = length(protein.names),
                        sets = rev(names(protein.names)),
                        keep.order = TRUE, mainbar.y.label=paste0("Intersection Size ", intersectionType)))
    dev.off()
  }, diaWorkflowResults.names=diaWorkflowResults.names)
}


plotUpsetPlot(diaWorkflowResults, selectedDilutions=c("Lymph nodes", "Lymph nodes + 1:25 E.coli", 
                                                          "Lymph nodes + 1:12 E.coli", "Lymph nodes + 1:6 E.coli"), str="allDilutions") 

plotUpsetPlot(diaWorkflowResults, selectedDilutions=c("Lymph nodes + 1:25 E.coli", 
                                                          "Lymph nodes + 1:12 E.coli"), str="12and25Only") 



########################################################################

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# TODO Adjust PCA Plot to having SkylineMSstats and SkylineDIANN

plotPCAGrid <- function(diaWorkflowResults, filename, labelled=FALSE) {
  pcaPlot.lst <- list()
  pcaPlot.lst <- lapply(1:length(diaWorkflowResults), function(idx){
    print(names(diaWorkflowResults)[idx])
    df <- diaWorkflowResults[[idx]]
    # Remove Sample 1-6_028 as it has too many missing values
    df <- df[, colnames(df) != "1-6_028"]
    t.df <- t(df)
    # For 1-6_028 in "Skyline_MaxQuant" are all protein itensities NA --> Remove all empty samples
    t.df <- t.df[rowSums(is.na(t.df)) != ncol(t.df),]
    
    # Remove proteins with only NAs
    t.df <- t.df[, colSums(is.na(t.df)) != nrow(t.df)]
    
    pca <- pcaMethods::pca(t.df, method="nipals", center = TRUE, maxSteps=5000)
    plot(pca)
    df2 <- merge(pcaMethods::scores(pca), t.df, by=0)
    df2 <- addDilutionColumn(df2, "Row.names")
    
    df2$Row.names <- as.numeric(gsub("1-6_|1-12_|1-25_", "", df2$Row.names))
    
    gg <- ggplot(df2, aes(PC1, PC2)) +
      xlab(paste("PC1 (", round(pca@R2[1], 3) * 100, "%)")) +
      ylab(paste("PC2 (", round(pca@R2[2], 3) * 100, "%)")) +
      # ggtitle(names(diaWorkflowResults)[idx]) +
      coord_fixed() +
      theme_minimal() +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            legend.title = element_blank(), legend.position="bottom",
            #legend.justification = "left", 
            legend.direction = "vertical") 
    
    if (labelled){
      gg <- gg +       
        geom_text(
          size = 2,
          label = df2$Row.names,
          nudge_x = 0.25, nudge_y = 0.25,
          check_overlap = F,
          aes(colour = Row.names)
        ) +
        viridis::scale_color_viridis(option = "mako", direction = -1)
    } else {
      gg <- gg + geom_point(aes(shape=dilution, color=dilution), alpha = 0.5)
    }
    gg
    
  })
  
  names(pcaPlot.lst) <- names(diaWorkflowResults)
  
  p2_legend <- get_legend(pcaPlot.lst[[5]]) 
  
  grobArranged <- gridExtra::arrangeGrob(pcaPlot.lst[["DIANN_PROSIT_EDIA_GPF"]] + theme(legend.position="none"), 
                              pcaPlot.lst[[ "DIANN_DIANN_AI_GPF"]] + theme(legend.position="none"), 
                              pcaPlot.lst[["DIANN_MaxQuant"]] + theme(legend.position="none"),
                              pcaPlot.lst[["DIANN_MSFragger"]] + theme(legend.position="none"), 
                              pcaPlot.lst[["DIANN_Predicted"]] + theme(legend.position="none"),
                              pcaPlot.lst[["Skyline_PROSIT_EDIA_GPF"]] + theme(legend.position="none"), 
                              pcaPlot.lst[["Skyline_DIANN_AI_GPF"]] + theme(legend.position="none"), 
                              pcaPlot.lst[["Skyline_MaxQuant"]] + theme(legend.position="none"),
                              pcaPlot.lst[["Skyline_MSFragger"]] + theme(legend.position="none"), 
                              patchwork::plot_spacer() + theme_minimal(), 
                              patchwork::plot_spacer() + theme_minimal(),
                              pcaPlot.lst[["OSW_DIANN_AI_GPF"]] + theme(legend.position="none"), 
                              pcaPlot.lst[["OSW_MaxQuant"]] + theme(legend.position="none"), 
                              pcaPlot.lst[["OSW_MSFragger"]] + theme(legend.position="none"), 
                              patchwork::plot_spacer() + theme_minimal(),
                              pcaPlot.lst[["Spectronaut_PROSIT_EDIA_GPF"]] + theme(legend.position="none"), 
                              pcaPlot.lst[["Spectronaut_DIANN_AI_GPF"]] + theme(legend.position="none"), 
                              pcaPlot.lst[["Spectronaut_MaxQuant"]] + theme(legend.position="none"),
                              pcaPlot.lst[["Spectronaut_MSFragger"]] + theme(legend.position="none"), 
                              pcaPlot.lst[["Spectronaut_Predicted"]] + theme(legend.position="none"), 
                              nrow=4)
  
  
  g <- grid.arrange(rbind(tableGrob(t(c("PROSIT EDIA GPF", "DIANN AI GPF", "MaxQuant", "MSFragger", "Predicted")), theme = ttheme_minimal(), rows = ""), 
                          cbind(tableGrob(c("DIANN", "Skyline", "OSW", "Spectronaut"), theme = ttheme_minimal()), 
                                grobArranged,  size = "last"), size = "last"), p2_legend, 
                    nrow=2, heights=c(10, 2))
  
  ggsave(file=filename, g, width=15, height=12)
}
# dev.off()
#plotPCAGrid(diaWorkflowResults, filename="pcaNIPALS_allDIAWorkflows.pdf") 
#plotPCAGrid(diaWorkflowResults, filename="pcaNIPALS_allDIAWorkflows_labelled2.pdf", labelled=TRUE) 


diaWorkflowResults.normalized <- lapply(diaWorkflowResults, function(df){
  df <- df[, colnames(df) != "1-6_028"]
  df.norm <- limma::normalizeQuantiles(df)
  df.norm
})
#plotPCAGrid(diaWorkflowResults.normalized, filename="pcaNIPALS_allDIAWorkflows_QNnormalized.pdf") 

# SUPPLEMENTARY FIGURE S1
plotPCAGrid(diaWorkflowResults.normalized, filename="FigS1_pcaNIPALS_QNnormalized_labelled.svg", labelled=TRUE) 



########################################################################

# GET NUMBER OF INTERSECT AND UNION PROTEINS FOR EACH BOOTSTRAP DATASET FOR EACH DIA WORKFLOW

combinedProteinNames <- readRDS("combinedProteinNames.rds")
intersectProteinNames <- readRDS("intersectProteinNames.rds")

nEcoli.comb <- length(combinedProteinNames[grepl("_ECOLI", combinedProteinNames)])
# # [1] 2125
nHuman.comb <- length(combinedProteinNames[grepl("_HUMAN", combinedProteinNames)])
# # [1] 11533
# 
nEcoli.intersect <- length(intersectProteinNames[grepl("_ECOLI", intersectProteinNames)])
# # [1] 740
nHuman.intersect <- length(intersectProteinNames[grepl("_HUMAN", intersectProteinNames)])
# # [1] 4512


# Get No. of proteins for each DIA workflow which were included for .intersect and .combined evaluation measures
nProteinsPerDiaWorkflow <- lapply(diaWorkflowResults, function(diaWorkflow){
  matchesCombined <- unlist(lapply(row.names(diaWorkflow), function(x) {
    length(intersect(unlist(base::strsplit(x, ";")), combinedProteinNames))>0
  }))
  combined <- diaWorkflow[matchesCombined,]
  nCombinedEcoli <- length(row.names(combined)[grepl("_ECOLI", row.names(combined))])
  nCombinedHuman <- length(row.names(combined)[grepl("_HUMAN", row.names(combined))])
  
  
  matchesIntersect <- unlist(lapply(row.names(diaWorkflow), function(x) {
    length(intersect(unlist(base::strsplit(x, ";")), intersectProteinNames))>0
  }))
  
  intersect <- diaWorkflow[matchesIntersect,]
  nIntersectEcoli <- length(row.names(intersect)[grepl("_ECOLI", row.names(intersect))])
  nIntersectHuman <- length(row.names(intersect)[grepl("_HUMAN", row.names(intersect))])
  
  
  list(nCombinedTotal = nrow(combined), nCombinedEcoli = nCombinedEcoli, nCombinedHuman = nCombinedHuman, 
       combinedEcoliRatio=round(nCombinedEcoli/nrow(combined), 4),
       nIntersectTotal = nrow(intersect), nIntersectEcoli=nIntersectEcoli, nIntersectHuman=nIntersectHuman,
       intersectEcoliRatio=round(nIntersectEcoli/nrow(intersect), 4))
})

df <- data.frame(matrix(unlist(nProteinsPerDiaWorkflow), nrow=length(nProteinsPerDiaWorkflow), byrow=TRUE))
row.names(df) <- names(nProteinsPerDiaWorkflow)
colnames(df) <- names(nProteinsPerDiaWorkflow[[1]])
df$dia <- row.names(df)


# df[df$dia == "DIANN_DIANN_AI",]$dia <- "DIANN_Predicted"
# df[df$dia == "Spectronaut_DirectDIA",]$dia <- "Spectronaut_Predicted"
# 
# df$dia <-gsub("_", " ", df$dia)
# df <- df[order(df$dia),]
df <- df[,c(ncol(df),1:(ncol(df)-1))]

# Part of this table is depicted in SUPPLEMENTARY FIGURE S4 AND S5
write.csv(df, file = "numberOfProteinsIntersectCombined.csv", row.names = FALSE)


#####################################################################################################################
#####################################################################################################################
# BOOTSTRAP

# READ IN FILES AND FORMAT RESULT DATA FRAME
theme_set(theme_minimal())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

listfiles <- list.files(path = "20210125_CSVs", full.names = TRUE, pattern = "*.csv") 

listfiles.entrynumbers <- lapply(X = listfiles, FUN = function(x) {
  length(count.fields(x, skip = 1))
})
names(listfiles.entrynumbers) <- listfiles
listfiles.entrynumbers <- t(data.frame(listfiles.entrynumbers))

df.combined <- ldply(listfiles, readr::read_csv)

factor.cols <- c("dia", "normalization", "sparsityReduction", "statTest")
df.combined[,factor.cols] <- lapply(df.combined[,factor.cols], factor) 

lapply(df.combined, class)

df.combined$dia <- as.character(df.combined$dia)

##########################################################################
# `%notin%` <- Negate(`%in%`)

df.combined$normalization <- factor(df.combined$normalization, levels = c("unnormalized","median", "QN", "TRQN"))
df.combined$sparsityReduction <- factor(df.combined$sparsityReduction, levels = c("NoSR", "SR66", "SR90"))
df.combined$statTest <- factor(df.combined$statTest, levels = c("ttestVarEqual", "ttestVarNotEqual", "GLMgamma", "limma", "Wilcoxon", "SAM", "ROTS"))
df.combined$portionEcoli <- df.combined$nEcoli/df.combined$nAllProteins
df.combined$nAllProteins.pre <- df.combined$nEcoli.pre + df.combined$nHuman.pre
df.combined$portionEcoliProteins.pre <- df.combined$nEcoli.pre/df.combined$nAllProteins.pre
df.combined$portionPC2PC1.woNAs <- df.combined$prctPC2.woNAs/df.combined$prctPC1.woNAs

df.combined[df.combined$dia == "DIANN_DIANN_AI",]$dia <- "DIANN_Predicted"
df.combined[df.combined$dia == "Spectronaut_DirectDIA",]$dia <- "Spectronaut_Predicted"
df.combined$dia <- gsub("_", " ", df.combined$dia)


# Sample uniqueness vs. performance of statistical tests

indices <- readRDS("indices.rds")
duplicate.df <- data.frame(do.call(rbind, lapply(indices, function (x) c(nUniqueSamples=length(unique(x)), nSamples=length(x)))))
duplicate.df$bootstrap.dataset <- row.names(duplicate.df)
duplicate.df$portionUniqueSamples <- duplicate.df$nUniqueSamples/duplicate.df$nSamples
duplicate.df <- duplicate.df[, c(3, 1, 2, 4)]

df.combined <- merge(df.combined, duplicate.df, by = "bootstrap.dataset")

write.csv(df.combined, "20220125_benchmark_analysis_results.csv", row.names = FALSE)
saveRDS(df.combined, "20220125_benchmark_analysis_results.rds")

settings <- c("dia", "normalization", "sparsityReduction", "statTest") 


df.combined <- df.combined[, c("bootstrap.dataset", "dia", "normalization", "sparsityReduction", "statTest", 
                               "groupSize", "medianSampleVariance", "medianProteinVariance", 
                               "percNATotal",  "kurtosis.wNAs", "skewness.wNAs", 
                               "var.groups.ratio.wNAs", "prctPC1.woNAs", "prctPC2.woNAs", "portionPC2PC1.woNAs", "p.pauc_0.9_correctFALSE.combined",
                               "sensAtpVal005.combined", "regpValsProp.combined", "p.pauc_0.9_correctFALSE.intersect",  
                               "sensAtpVal005.intersect", "regpValsProp.intersect", "p.pauc_0.9_correctFALSE.diaWorkflow", 
                               "sensAtpVal005.diaWorkflow", "regpValsProp.diaWorkflow", 
                               "RMSEEcoli.intersect", "RMSEHuman.intersect", "RMSEHumanAndEcoli.intersect",
                               "RMSEEcoli.diaWorkflow", "RMSEHuman.diaWorkflow", "RMSEHumanAndEcoli.diaWorkflow",
                               "portionUniqueSamples")]

eval.measure <- "p.pauc_0.9_correctFALSE.diaWorkflow"

########################################################################

maxpAUC <- max(df.combined$p.pauc_0.9_correctFALSE.diaWorkflow)

selectedGroupSizes <- c(3, 6, 13, 23)
for (selectedGroupSize in selectedGroupSizes){
  for (selectedSR in unique(df.combined$sparsityReduction)){
    
    df.combined1 <- df.combined[df.combined$sparsityReduction == selectedSR & df.combined$groupSize==selectedGroupSize,]
    
    plots <- df.combined1 %>% dplyr::group_by(df.combined1$dia) %>%
      do(
        plots = ggplot(data = .) + aes(x=portionUniqueSamples, y=p.pauc_0.9_correctFALSE.diaWorkflow, color=statTest) +
          geom_point(alpha=0.1, size = 0.25) + ggtitle(paste(.$dia, selectedSR, selectedGroupSize, sep = ", ")) +
          guides(colour = guide_legend(override.aes = list(alpha = 1))) +
          theme_minimal() +
          geom_smooth(method=lm, aes(fill=statTest)) +
          ylim(0, maxpAUC)
      )
    
    pdf(file=paste0("S19_correlation_portionUniqueSamples_pAUC_", selectedSR, "only_groupSize", selectedGroupSize, ".pdf"))
    for (plot in plots$plots){
      print(plot)
    }
    dev.off()
  }
}
########################################################################
# library("scattermore")
# plot(df.combined$RMSEEcoli.diaWorkflow, df.combined$RMSEHuman.diaWorkflow)
# 
# df.combined.LibSoftware <- addSoftwareLibraryColumns(df.combined, sepSpace=TRUE)
# 
# print(ggplot(df.combined.LibSoftware, aes(x=RMSEHumanAndEcoli.diaWorkflow, y=p.pauc_0.9_correctFALSE.diaWorkflow, col=diaSoftware)) +
#         geom_scattermore(alpha=0.01)+
#         theme(legend.title = element_blank()) +
#         guides(colour = guide_legend(override.aes = list(alpha = 1))))
# 
# 
# 
# 
# pdf(file=paste0("S_correlation_RMSE_pAUC90_2.pdf"))
# print(ggplot(df.combined.LibSoftware, aes(x=RMSEHumanAndEcoli.diaWorkflow, y=p.pauc_0.9_correctFALSE.diaWorkflow, col=sparsityReduction)) +
#         geom_point(alpha=0.05, pch='.') + 
#         theme(legend.title = element_blank()) +
#         facet_grid(diaSoftware ~ diaLibrary, scales = "fixed") +
#         guides(colour = guide_legend(override.aes = list(alpha = 1))))
# dev.off()
# 
# for (setting in setdiff(settings, "dia")) {
#   pdf(file=paste0("S_correlation_RMSE_pAUC90_", setting, ".pdf"), height = 10, width = 10)
#   print(ggplot(df.combined.LibSoftware, aes(x=RMSEHumanAndEcoli.diaWorkflow, y=p.pauc_0.9_correctFALSE.diaWorkflow, col=.data[[setting]])) +
#           geom_point(alpha=0.02, pch='.') + 
#           theme(legend.title = element_blank()) +
#           facet_grid(diaSoftware ~ diaLibrary, scales = "fixed") +
#           guides(colour = guide_legend(override.aes = list(alpha = 1))))
#   dev.off() 
# }

########################################################################
# Performance without separation by settings

plotViolinPlotsCombSetting <- function(setting, df.combined, eval.measure) {
  ggplot.violin <- ggplot(df.combined, aes(y = forcats::fct_rev(get(setting)), x = get(eval.measure))) + 
    ylab("") +
    xlab(eval.measure) +
    geom_violin(scale="count") +
    geom_boxplot(width=0.08, outlier.size=0.5, outlier.alpha=0.4) +
    geom_vline(aes(xintercept = median(get(eval.measure), na.rm = TRUE)), colour = 'red') 
  
  ggplot.violin
}

plotViolinPlotsCombSettingHorizontal <- function(setting, df.combined, eval.measure) {
  ggplot.violin <- ggplot(df.combined, aes(x = get(setting), y = get(eval.measure))) + 
    xlab("") +
    ylab(eval.measure) +
    geom_violin(scale="count") +
    geom_boxplot(width=0.08, outlier.size=0.5, outlier.alpha=0.4) +
    geom_hline(aes(yintercept = median(get(eval.measure), na.rm = TRUE)), colour = 'red') +
    theme(strip.placement = "outside", axis.text.x = element_text(angle = 55, hjust = 1))
  
  ggplot.violin
}

# FIGURE 4
eval.measures <- c("p.pauc_0.9_correctFALSE.diaWorkflow", "p.pauc_0.9_correctFALSE.intersect", "p.pauc_0.9_correctFALSE.combined")

for (eval.measure in eval.measures){
  plotsSettings <- lapply(settings, plotViolinPlotsCombSettingHorizontal, df.combined=df.combined, eval.measure=eval.measure)
  pdf(file = paste0("violinPlots_settingsCombined_horizontal_", eval.measure, ".pdf"), height = 5, width = 18)
  print(cowplot::plot_grid(plotsSettings[[1]], plotsSettings[[2]], plotsSettings[[4]], plotsSettings[[3]], align = "h", nrow=1))
  dev.off()
}

for (eval.measure in eval.measures){
  plotsSettings <- lapply(settings, plotViolinPlotsCombSetting, df.combined=df.combined, eval.measure=eval.measure)
  pdf(file = paste0("Fig4CAndD_violinPlots_settingsCombined_", eval.measure, ".pdf"), height = 22, width = 6)
  print(cowplot::plot_grid(plotsSettings[[1]], plotsSettings[[2]], plotsSettings[[4]], plotsSettings[[3]], align = "v", ncol=1))
  dev.off()
}

eval.measure <- "p.pauc_0.9_correctFALSE.diaWorkflow"
selectedSR <- "NoSR"

df.combined.SelectedSRonly <- df.combined[df.combined$sparsityReduction==selectedSR,]
plotsSettings.SelectedSR <- lapply(setdiff(settings, "sparsityReduction"), plotViolinPlotsCombSetting,df.combined=df.combined.SelectedSRonly, eval.measure=eval.measure)

svg(file = paste0("violinPlots_settingsCombined_", selectedSR, "_", eval.measure,".svg"), height = 18, width = 5)
cowplot::plot_grid(plotlist=plotsSettings.SelectedSR, align = "v", ncol=1)
dev.off()

eval.measure <- "p.pauc_0.9_correctFALSE.diaWorkflow"
svg(file = paste0("Fig4B_violinPlots_settingsCombined_", eval.measure, "_SRs.svg"), height = 5, width = 5)
ggplot(df.combined, aes(x = sparsityReduction, y = get(eval.measure))) + 
  xlab("") +
  ylab(eval.measure) +
  geom_violin(scale="count") +
  geom_boxplot(width=0.08, outlier.size=0.5, outlier.alpha=0.4) +
  geom_hline(aes(yintercept = median(get(eval.measure), na.rm = TRUE)), colour = 'red') 
dev.off()

########################################################################

df.combined[,"groupSize"] <- as.factor(df.combined[,"groupSize"])
settChars <- c(settings[-1], "groupSize")

getViolinGrid <- function(settChar, eval.measure, df, str=""){
  ggplot.violin <- ggplot(df, aes(x = get(settChar), y = get(eval.measure))) + 
    xlab(settChar) + ylab(eval.measure) +
    geom_violin(scale="count") +
    geom_boxplot(width=0.1, outlier.size=0.5, outlier.alpha=0.4) +
    facet_grid(diaSoftware ~ diaLibrary, scales = "fixed") +
    theme_bw() +
    xlab("") +
    theme(strip.placement = "outside", axis.text.x = element_text(angle = 55, hjust = 1)) +
    geom_hline(aes(yintercept = median(get(eval.measure), na.rm = TRUE)), colour = 'red') 
  ggplot2::ggsave(file = paste0(str, eval.measure, "_", settChar,".svg"), plot = ggplot.violin, device = "svg", height = 8, width = 14)
}

df.combined.LibSoftware <- df.combined
df.combined.LibSoftware <- addSoftwareLibraryColumns(df.combined.LibSoftware, sepSpace=TRUE)

df.combined.LibSoftware.NoSROnly <- df.combined.LibSoftware[df.combined.LibSoftware$sparsityReduction=="NoSR",]

# SUPPLEMENTARY FIGURE S6
getViolinGrid("normalization","p.pauc_0.9_correctFALSE.diaWorkflow", df.combined.LibSoftware.NoSROnly, str="FigS13_normalization_violin_NoSROnly_")

# SUPPLEMENTARY FIGURE S7
getViolinGrid("statTest","p.pauc_0.9_correctFALSE.diaWorkflow", df.combined.LibSoftware.NoSROnly, str="FigS14_statTest_violin_NoSROnly_")

df.combined.LibSoftware.NoSROnly$groupSize <- as.numeric(as.character(df.combined.LibSoftware.NoSROnly$groupSize))

# SUPPLEMENTARY FIGURE S17
for (eval.measure2 in c("p.pauc_0.9_correctFALSE.diaWorkflow", 
                        "p.pauc_0.9_correctFALSE.intersect", 
                        "p.pauc_0.9_correctFALSE.combined",
                        "sensAtpVal005.diaWorkflow", 
                        "sensAtpVal005.intersect", 
                        "sensAtpVal005.combined" 
)){
  ggplot.medianLine.groupSize.statTest <- ggplot(df.combined.LibSoftware.NoSROnly, aes(x = groupSize, y = get(eval.measure2))) +
    xlab("groupSize") + ylab(eval.measure2) +
    stat_summary(fun="median", geom="line", aes(color=statTest)) +
    facet_grid(diaSoftware ~ diaLibrary, scales = "fixed") +
    theme_bw() +
    xlab("") +
    geom_hline(aes(yintercept = median(get(eval.measure2), na.rm = TRUE)), colour = 'red') +
    theme(legend.title = element_blank(), legend.position="bottom",
          legend.justification = "left")
  
  if (grepl("sensAtpVal005.", eval.measure2)){
    ggplot.medianLine.groupSize.statTest <- ggplot.medianLine.groupSize.statTest + geom_hline(aes(yintercept = 0.8), colour = 'blue')
  }
  ggsave(paste0("FigS17_ggplot_medianLine_groupSize_NoSROnly_", eval.measure2, ".svg"), ggplot.medianLine.groupSize.statTest, height = 8, width = 10)
}

df.combined.LibSoftware.selected <- df.combined.LibSoftware[, c("dia",
                                                                "sparsityReduction",
                                                                "groupSize",
                                                                "statTest",
                                                                "p.pauc_0.9_correctFALSE.combined",
                                                                "p.pauc_0.9_correctFALSE.intersect",
                                                                "p.pauc_0.9_correctFALSE.diaWorkflow")]

# selectedDia <- "DIANN DIANN AI GPF"
for (selectedDia in unique(df.combined.LibSoftware.selected$dia)){
  df.combined.LibSoftware.selected.dia <- df.combined.LibSoftware.selected[df.combined.LibSoftware.selected$dia == selectedDia,]
  
  df.combined.LibSoftware.selected.dia$groupSize <- factor(df.combined.LibSoftware.selected.dia$groupSize)
  df.combined.LibSoftware.selected.dia2 <- reshape2::melt(df.combined.LibSoftware.selected.dia)
  df.combined.LibSoftware.selected.dia2$proteinList <- sub('.*\\.', '', df.combined.LibSoftware.selected.dia2$variable)
  df.combined.LibSoftware.selected.dia2$groupSize <- as.numeric(as.character(df.combined.LibSoftware.selected.dia2$groupSize))
  
  
  ggplot.medianLine.groupSize.SRProteinSetstatTest <- ggplot(df.combined.LibSoftware.selected.dia2, aes(x = groupSize, y = value)) +
    xlab("Group size") + ylab("pAUC") +
    stat_summary(fun="median", geom="line", aes(color=statTest)) +
    facet_grid(proteinList ~ sparsityReduction, scales = "fixed") +
    theme_bw() +
    geom_hline(aes(yintercept = median(value, na.rm = TRUE)), colour = 'red') +
    theme(legend.title = element_blank(), legend.position="bottom",
          legend.justification = "left")
  
  ggsave(paste0("medianLine_SRProteinSetstatTest_", gsub(" ", "_", selectedDia), ".svg"), ggplot.medianLine.groupSize.SRProteinSetstatTest, height = 8, width = 10)
}

########################################################################

data.characteristics.sel <- c("groupSize",  
                              "medianSampleVariance", "medianProteinVariance", 
                              "percNATotal", "kurtosis.wNAs", "skewness.wNAs", 
                              "var.groups.ratio.wNAs", "prctPC1.woNAs", "prctPC2.woNAs")

eval.measures.sel <- sort(c("p.pauc_0.9_correctFALSE.combined", 
                            "sensAtpVal005.combined", "regpValsProp.combined", 
                            "p.pauc_0.9_correctFALSE.intersect",
                            "sensAtpVal005.intersect", "regpValsProp.intersect", 
                            "p.pauc_0.9_correctFALSE.diaWorkflow", 
                            "sensAtpVal005.diaWorkflow", "regpValsProp.diaWorkflow", 
                            "RMSEEcoli.diaWorkflow", "RMSEHuman.diaWorkflow", 
                            "RMSEEcoli.intersect", "RMSEHuman.intersect"))

sett <- "dia"  
sel <- c(data.characteristics.sel, eval.measures.sel)

df.combined.pairsplotsel <- df.combined[df.combined$sparsityReduction=="NoSR", c(sett, sel)]
df.combined.pairsplotsel$groupSize <- as.numeric(as.character(df.combined.pairsplotsel$groupSize))

# SUPPLEMENTARY FIGURE S16
for (dia in sort(unique(df.combined.pairsplotsel$dia))){
  df.combined.pairsplotsel.dia <- df.combined.pairsplotsel[df.combined.pairsplotsel$dia==dia,]
  df.combined.pairsplotsel.dia$dia <- NULL
  
  M <- cor(df.combined.pairsplotsel.dia, use="pairwise.complete.obs") # "pairwise.complete.obs" as there are NAs in regVals for some dias
  
  svg(file = paste0("FigS16_corrplot_charact_NoSR_", dia, ".svg"), width = 15, height = 15)
  corrplot::corrplot.mixed(M,  upper = "ellipse", lower = "number",
                 tl.pos = "lt", tl.col = "black", main=dia, mar=c(0,0,2,0))
  dev.off()
}

df.combined.pairsplotsel.long <- reshape2::melt(df.combined.pairsplotsel[, c("dia", setdiff(data.characteristics.sel, "groupSize"))])
df.combined.pairsplotsel.long <- addSoftwareLibraryColumns(df.combined.pairsplotsel.long, sepSpace = TRUE)

# SUPPLEMENTARY FIGURE S15
ggplot.charact <- ggplot(df.combined.pairsplotsel.long, aes(forcats::fct_rev(dia), value)) +
  geom_boxplot(aes(fill=diaSoftware), alpha=0.5, outlier.size=0.5) +
  coord_flip() +
  xlab("") +
  ylab("") +
  facet_wrap( ~ variable, scales = "free", ncol=2) +
  ggplot2::theme_bw() +
  theme(legend.title = element_blank(), axis.text.y = element_text(hjust=0), legend.position="bottom",
        legend.justification = "left", legend.direction = "vertical") +
  guides(fill = guide_legend(reverse = TRUE))
ggsave(file=paste0("FigS15_charact_bootstrapDatasets_NoSROnly.svg"), ggplot.charact, height=14, width=11)

########################################################################

sett <- "dia"
sel <- c("bootstrap.dataset", "groupSize", settings)

eval.measures <- c("p.pauc_0.9_correctFALSE.diaWorkflow", "p.pauc_0.9_correctFALSE.intersect", "p.pauc_0.9_correctFALSE.combined")
for (eval.measure in eval.measures){
  df.combined.eval <- df.combined[, c(sel, eval.measure)]
  
  df.combined.eval$concatenated <- do.call(paste, c(df.combined.eval[, setdiff(sel, sett)], sep = "_"))
  df.combined.eval <- df.combined.eval[, c(sett, "concatenated", eval.measure)]
  
  df.combined.eval.wide <- tidyr::spread(df.combined.eval, key = concatenated, value = get(eval.measure))
  row.names(df.combined.eval.wide) <- df.combined.eval.wide[,sett]
  df.combined.eval.wide[,sett] <- NULL
  df.combined.eval.wide <- t(df.combined.eval.wide)
  
  M <- cor(df.combined.eval.wide, use="pairwise.complete.obs")
  
  # SUPPLEMENTARY FIGURE S8
  svg(file = paste0("FigS8_corrplot_", sett, "_",eval.measure,".svg"), width = 15, height = 15)
  corrplot.mixed(M,  upper = "ellipse", lower = "number",
                 tl.pos = "lt", tl.col = "black")
  dev.off()
  
  svg(file = paste0("corrplot2_", sett, "_",eval.measure,".svg"), width = 15, height = 15)
  corrplot(M, method = "number", tl.col="black")
  dev.off()
  
  svg(file = paste0("networkplot_", sett, "_",eval.measure,".svg"), width = 20, height = 20)
  print(df.combined.eval.wide %>% 
          corrr::correlate() %>% 
          corrr::network_plot())
  dev.off()
}

########################################################################
# Ranking

df.combined2 <- df.combined
df.combined2 <- addSoftwareLibraryColumns(df.combined2, sepSpace = TRUE)

diaSpecifications <- c("dia", "diaSoftware", "diaLibrary")
settings <- c("statTest", "sparsityReduction", "normalization")
proteinLists <- c("diaWorkflow", "intersect", "combined")

for (setting in settings){
  for (diaSpecification in diaSpecifications){
    bump.lst <- list()
    for (proteinList in proteinLists) {
      ranks.df <- df.combined2 %>% 
        dplyr::group_by(bootstrap.dataset, !!as.name(setting), dia, diaSoftware, diaLibrary) %>% 
        dplyr::summarise(mean.pAUC = mean(!!as.name(paste0("p.pauc_0.9_correctFALSE.", proteinList)))) %>% 
        dplyr::group_by(bootstrap.dataset, dia) %>%
        dplyr::mutate(ranking = rank(- mean.pAUC)) %>% 
        dplyr::group_by(!!as.name(setting), !!as.name(diaSpecification)) %>% 
        dplyr::summarise(mean.ranking = mean(ranking))
      
      bump.lst <- list.append(bump.lst, ranks.df)
      
      if (diaSpecification == "dia"){
        width = 8
      } else if (diaSpecification == "diaLibrary") {
        width = 7
      } else {
        width = 7
      }
      
      bump.chart <- ggplot(data = ranks.df, aes(x = get(setting), y = mean.ranking, group = get(diaSpecification), alpha=0.5)) +
        geom_line(aes(color = get(diaSpecification)), size = 1) +
        geom_point(aes(color = get(diaSpecification)), size = 2) +
        scale_y_reverse(limits=c(length(unique(ranks.df[[setting]],0)), 1), breaks = 1:nrow(ranks.df)) +
        ggtitle(paste0(diaSpecification, ", ", setting, ", ", proteinList)) +
        theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45)) +
        # xlab(setting) +
        xlab("") +
        scale_alpha(guide = 'none') +
        guides(colour = guide_legend(override.aes = list(alpha = 0.5)))
      ggsave(paste0("S18_bumpchart_", diaSpecification, "_", setting, "_", proteinList, ".svg"), bump.chart, width=width, height=5)
      
    }
    
    meanOverProteinLists.df <- data.frame(bump.lst[[1]][,1:2], mean.mean.ranking=rowMeans(data.frame(bump.lst[[1]][,"mean.ranking"], bump.lst[[2]][,"mean.ranking"], bump.lst[[3]][,"mean.ranking"])))
    bump.chart <- ggplot(data = meanOverProteinLists.df, aes(x = get(setting), y = mean.mean.ranking, group = get(diaSpecification), alpha=0.5)) +
      geom_line(aes(color = get(diaSpecification)), size = 1) +
      geom_point(aes(color = get(diaSpecification)), size = 2) +
      scale_y_reverse(limits=c(length(unique(ranks.df[[setting]],0)), 1), breaks = 1:nrow(meanOverProteinLists.df)) +
      ggtitle(paste0(diaSpecification, ", ", setting)) +
      theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45)) +
      # xlab(setting) +
      xlab("") +
      scale_alpha(guide = 'none') +
      guides(colour = guide_legend(override.aes = list(alpha = 0.5)))
    ggsave(paste0("S18_bumpchart_", diaSpecification, "_", setting, "_meanOverproteinLists.svg"), bump.chart, width=width, height=5)
  }
}

session <- sessionInfo()
sink("sessionInfo_visualizations.txt")
print(session)
sink()