library(data.table)
library(gridExtra)
library(ggplot2)
theme_set(theme_bw())
library(RColorBrewer)
library(ggridges)
library(rlist)
library(tidyr)
library(plyr)
library(dplyr)
library(UpSetR)
library(viridis)
library(ggplot2)
library(patchwork)
library(grid)
library(readr) # Dev version: devtools::install_github("tidyverse/readr")
library(ggplot2)
library(corrplot)
library(RColorBrewer)
library(reshape2)
library(GGally)
library(corrplot)
library(corrr)
library(forcats)
library(pcaMethods)
library(limma)
library(cowplot)

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
  
  # Replace DIANN_DIANN_AI by DIANN_Predicted and Spectronaut_DirectDIA by Spectronaut_Predicted
  if ("DIANN_DIANN_AI" %in% unique(df.combined.LibSoftware$dia)){
    df.combined.LibSoftware[df.combined.LibSoftware$dia == "DIANN_DIANN_AI",]$dia <- "DIANN_Predicted"
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

#diaWorkflowResults <- lapply("Converted_to_same_format.RData", function(x) mget(load(x)))
diaWorkflowResults <- lapply("DIAsoftwareOutputProteinLevel.RData", function(x) mget(load(x)))
diaWorkflowResults <- unlist(diaWorkflowResults,recursive=FALSE)
diaWorkflowResults <- list.remove(diaWorkflowResults, ".Random.seed")

formatMSTable <- function(x) {
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
}

addLeadingZeros <- function(el){
  if (length(el) == 1){
    elFinal <- formatC(as.numeric(el),width=3,format='f',digits=0,flag='0')
  } else if (length(el) == 2){
    elZerosAdded <- formatC(as.numeric(el[2]),width=3,format='f',digits=0,flag='0')
    elFinal <- paste0(el[1], "_", elZerosAdded)
  }
  # print(elFinal)
  elFinal
}

generateWideFormattedDf <- function(df, log=TRUE) {
  # remove duplicates from multiple precursors carrying the same quantitative information for a certain protein
  df <- dplyr::distinct(df)
  
  df <- reshape(df, idvar="Protein.Names", timevar="Sample",  direction="wide")
  colnames(df) <- sub("PG.Quantity.", "", colnames(df))
  
  row.names(df) <- df$Protein.Names
  df$Protein.Names <- NULL
  df[df == 0] <- NA
  
  if (log) df <- log2(df)
  
  return(df)
}

cleanData <- function(x, log = TRUE) {
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
  x <- generateWideFormattedDf(x, log)
  
  # Remove empty rows
  x <- x[rowSums(is.na(x)) != ncol(x),]
  x
}

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

# UNIFY OUTPUT TABLES
diaWorkflowResults <- lapply(diaWorkflowResults, formatMSTable)

# Change Protein names in OSW_MaxQuant from format Entry to Entry Name
OSW_MaxQuant_unique <- unique(diaWorkflowResults[["OSW_MaxQuant"]]$Protein.Names)

# Derived translation table from https://www.uniprot.org/uploadlists/
mappedProteins <- read.table("OSW_MaxQuant_unique_mappedProteins_20210204.tsv", sep = "\t", header = TRUE)
diaWorkflowResults[["OSW_MaxQuant"]]$Protein.Names <- mappedProteins$Entry.name[match(diaWorkflowResults[["OSW_MaxQuant"]]$Protein.Names, mappedProteins$Entry)]

# Only keep Entry Name from Protein Name in Skyline_MaxQuant, Skyline_PROSIT_EDIA_GPF, Skyline_DIANN_AI_GPF
diaWorkflowResults[["Skyline_MaxQuant"]]$Protein.Names <- sub(".*\\|", "", diaWorkflowResults[["Skyline_MaxQuant"]]$Protein.Names)
diaWorkflowResults[["Skyline_PROSIT_EDIA_GPF"]]$Protein.Names <- sub(".*\\|", "", diaWorkflowResults[["Skyline_PROSIT_EDIA_GPF"]]$Protein.Names)
diaWorkflowResults[["Skyline_DIANN_AI_GPF"]]$Protein.Names <- sub(".*\\|", "", diaWorkflowResults[["Skyline_DIANN_AI_GPF"]]$Protein.Names)
diaWorkflowResults[["Spectronaut_MSFragger"]]$Protein.Names <- gsub("'", '', sub(".*\\|", "", diaWorkflowResults[["Spectronaut_MSFragger"]]$Protein.Names))
diaWorkflowResults[["Skyline_MSFragger"]]$Protein.Names <- sub(".*\\|", "", diaWorkflowResults[["Skyline_MSFragger"]]$Protein.Names)

diaWorkflowResults[["OSW_MaxQuant"]]$Sample <- unlist(lapply( strsplit(unique(diaWorkflowResults[["OSW_MaxQuant"]]$Sample),split='_', fixed=TRUE), addLeadingZeros))
diaWorkflowResults[["OSW_DIANN_AI_GPF"]]$Sample <- unlist(lapply( strsplit(unique(diaWorkflowResults[["OSW_DIANN_AI_GPF"]]$Sample),split='_', fixed=TRUE), addLeadingZeros))
diaWorkflowResults[["OSW_MSFragger"]]$Sample <- unlist(lapply( strsplit(unique(diaWorkflowResults[["OSW_MSFragger"]]$Sample),split='_', fixed=TRUE), addLeadingZeros))

diaWorkflowResults.log <- lapply(diaWorkflowResults, cleanData, log=TRUE)
diaWorkflowResults.log <- diaWorkflowResults.log[order(names(diaWorkflowResults.log))]

design <- getDesign(diaWorkflowResults.log[[1]])

# Bring column names in same order as design
diaWorkflowResults.log <- lapply(diaWorkflowResults.log, function(x, design){
  x <- x[row.names(design)]
  x
}, design=design)

names(diaWorkflowResults.log)[names(diaWorkflowResults.log) == "DIANN_DIANN_AI"] <- "DIANN_Predicted"
names(diaWorkflowResults.log)[names(diaWorkflowResults.log) == "Spectronaut_DirectDIA"] <- "Spectronaut_Predicted"
names(diaWorkflowResults.log) <-  gsub("_", " ", names(diaWorkflowResults.log))
diaWorkflowResults.log <- diaWorkflowResults.log[order(names(diaWorkflowResults.log))]

########################################################################

stacked <- lapply(1:length(diaWorkflowResults.log), function(idx){
  x <- diaWorkflowResults.log[[idx]]
  x$ProteinName <- row.names(x)
  x <- reshape2::melt(x)
  x$Species <- "Human"
  x[grepl("_ECOLI", x$ProteinName),]$Species <- "E. coli"
  x <- addDilutionColumn(x, "variable") 
  x$dia <- names(diaWorkflowResults.log)[idx]
  colnames(x) <- c("proteinName", "sample", "intensity", "species", "dilution", "dia")
  x
})

stacked.df <- do.call(rbind, stacked)
stacked.df$dia <- gsub("_", " ", stacked.df$dia )
stacked.df <- addSoftwareLibraryColumns(stacked.df, sepSpace = TRUE)
stacked.df$dilution <- factor(stacked.df$dilution,
                              levels = c("Lymph nodes", "Lymph nodes + 1:25 E.coli", "Lymph nodes + 1:12 E.coli", "Lymph nodes + 1:6 E.coli"))

stacked.nonNA <- stacked.df %>% filter(!is.na(intensity)) %>% dplyr::group_by(sample, dilution, dia) %>% 
  dplyr::summarise(NumberProteins = n())
stacked.nonNA.meanMedians <- stacked.nonNA %>%  dplyr::group_by(dilution, dia) %>% 
  dplyr::summarise(median = median(NumberProteins), mean = mean(NumberProteins))

# FIGURE 2
ggplot.intensity.violinBoxplot <- ggplot(stacked.nonNA,aes(forcats::fct_rev(dia), NumberProteins, fill=forcats::fct_rev(dilution))) +
  geom_boxplot(alpha=0.5, outlier.size=0.5) +
  coord_flip() +
  xlab("") +
  ylab("# Proteins") +
  ggplot2::theme_minimal() +
  ggplot2::scale_fill_brewer(palette="Greens", direction=-1) +
  theme(legend.title = element_blank(), axis.text.y = element_text(hjust=0), legend.position="bottom",
        legend.justification = "left", legend.direction = "vertical") +
  guides(fill = guide_legend(reverse = TRUE)) 
ggsave(file="Fig2_ProteinNumbers.pdf", ggplot.intensity.violinBoxplot, height=8, width=6)

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
ggsave(file="Fig2_intensityDistributions.pdf", ggplot.intensityFacet, height=8, width=6)

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
  group_by(proteinName, dilution, dia, diaSoftware, diaLibrary, species) %>% 
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
ggsave(file="ggplot.intensityvariance.boxplot_Ecoli.pdf", ggplot.var, height=8, width=6)

# FIGURE 2
ggplot.var.mins15 <- ggplot.var + xlim(-15, NA) 
ggsave(file="Fig2_proteinVariance_Ecoli_xlimMinus15.pdf", ggplot.var.mins15, height=8, width=6)


########################################################################

proteinMeansNAsPlot.lst <- lapply(1:length(diaWorkflowResults.log), function(idx){
  dia <- names(diaWorkflowResults.log)[idx]
  print(dia)
  df <- diaWorkflowResults.log[[idx]]
  
  proteinMeans <- rowMeans(df, na.rm=TRUE)
  percNAsInProtein <- (rowSums(is.na(df))/ncol(df)) * 100
  protein.df <- data.frame(percNAsInProtein, proteinMeans, dia, proteinName=row.names(df))
  protein.df
})

getSampleNADf <-  function(idx, diaWorkflowResults.log, species=""){
  dia <- names(diaWorkflowResults.log)[idx]
  print(dia)
  df <- diaWorkflowResults.log[[idx]]
  
  if (species != ""){
    df <- df[grepl(paste0("_", species), row.names(df)),]
  }
  
  sampleMeans <- colMeans(df, na.rm=TRUE)
  percNAsInSample <- (colSums(is.na(df))/nrow(df)) * 100
  sample.df <- data.frame(percNAsInSample, sampleMeans, dia, sampleName=colnames(df))
  
  sample.df
}

getSampleNADfPerSpecies <- function(diaWorkflowResults.log, species="") {
  sampleMeansNAsPlot.humanAndEcoli.df <- do.call(rbind, 
                                                 lapply(1:length(diaWorkflowResults.log), 
                                                        getSampleNADf, 
                                                        diaWorkflowResults.log=diaWorkflowResults.log, 
                                                        species=species))
  sampleMeansNAsPlot.humanAndEcoli.df <- addDilutionColumn(sampleMeansNAsPlot.humanAndEcoli.df, "sampleName") 
  sampleMeansNAsPlot.humanAndEcoli.df
}

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
  pdf(file=paste0("Fig3A_proteinMeansNACorrelation_", gsub(" |\\.", "_", species),".pdf"), width=8, height=6)
  if (species == "E. coli & Human"){
    
    print(ggplot(df, aes(x = percNAsInProtein, y = proteinMeans, color=forcats::fct_rev(species))) + 
            geom_point(size=0.5, alpha=.05) +
            facet_grid(diaSoftware ~ diaLibrary, scales = "free") +
            scale_color_manual(values=c(brewer.pal(11, "RdYlBu")[10], brewer.pal(11, "RdYlBu")[2])) +
            xlab("% Missing values") +
            ylab("Protein mean") +
            ggtitle(species) +
            theme(legend.title = element_blank(), axis.text.y = element_text(hjust=0), legend.position="bottom",
                  legend.justification = "left", legend.direction = "vertical") +
            guides(colour = guide_legend(override.aes = list(alpha = 1))))
    
  } else {
    print(ggplot(df, aes(x = percNAsInProtein, y = proteinMeans)) + 
            geom_point(size=0.5, alpha=.05) +
            facet_grid(diaSoftware ~ diaLibrary, scales = "free") +
            xlab("% Missing values") +
            ylab("Protein mean") +
            ggtitle(species) +
            theme(legend.title = element_blank(), axis.text.y = element_text(hjust=0), legend.position="bottom",
                  legend.justification = "left", legend.direction = "vertical") )
  }
  dev.off()
  
  if (species == "E. coli & Human"){
    df2 <- getSampleNADfPerSpecies(diaWorkflowResults.log, species="")
    df2 <- addSoftwareLibraryColumns(df2, sepSpace = TRUE) 
    
    # FIGURE 3B
    pdf(file=paste0("Fig3B_sampleMeansNACorrelation_", gsub(" |\\.", "_", species),".pdf"), width=8, height=6)
    print(ggplot(df2, aes(x = percNAsInSample, y = sampleMeans, color=forcats::fct_rev(dilution))) + 
            geom_point(size=0.9, alpha=0.5, stroke = 0.1) +
            facet_grid(diaSoftware ~ diaLibrary,  scales = "free") +
            xlab("% Missing values") +
            ylab("Sample mean") +
            ggtitle(species) +
            theme(legend.title = element_blank(), axis.text.y = element_text(hjust=0), legend.position="bottom",
                  legend.justification = "left", legend.direction = "vertical") +
            guides(color = guide_legend(override.aes = list(alpha = 1))) )
    dev.off()
  }
}

########################################################################

# design <- getDesign(diaWorkflowResults.log[[1]])
# 
# # Bring column names to same order as design
# diaWorkflowResults.log <- lapply(diaWorkflowResults.log, function(x, design){
#   x <- x[row.names(design)]
#   x
# }, design=design)

diaWorkflowResults.log.rbind <- do.call(rbind, diaWorkflowResults.log)

mtxTwoGroups <- diaWorkflowResults.log.rbind[,grepl("1-25|1-12", colnames(diaWorkflowResults.log.rbind))]
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

# SUPPLEMENTARY FIGURE S1
pdf(file="FigS1_FC_ecoli.pdf", width=12, height=8)
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
    print(dilution)
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
  
  if (intensityCombOrMean == "intensityProteinMean"){
    df$ProteinNames <- row.names(df)
    df <-reshape2::melt(df, id.vars = c("ProteinNames","dia"))
  }
  df
}

getCombinedIntensitiesPerDIAWorkflow <- function(diaWorkflowResults.log, humanOrEcoli=""){
  intensitiesPerDilution <- lapply(diaWorkflowResults.log, getIntensitiesPerDilution, humanOrEcoli=humanOrEcoli)
  intensitiesPerDIAWorkflow <- lapply(1:length(intensitiesPerDilution), getCombIntensitiesMeans, 
                                      intensitiesPerDilution=intensitiesPerDilution,intensityCombOrMean= "intensityCombined")
  intensitiesPerDIAWorkflow.df <- do.call("rbind", intensitiesPerDIAWorkflow)
  colnames(intensitiesPerDIAWorkflow.df) <- c("intensity", "dilution", "dia")
  intensitiesPerDIAWorkflow.df <- addSoftwareLibraryColumns(intensitiesPerDIAWorkflow.df, sepSpace = TRUE)
  intensitiesPerDIAWorkflow.df
}

getProteinMeansIntensitiesPerDIAWorkflow <- function(diaWorkflowResults.log, humanOrEcoli=""){
  proteinMeansPerDilution <- lapply(diaWorkflowResults.log, getIntensitiesPerDilution, humanOrEcoli=humanOrEcoli)
  proteinMeansPerDIAWorkflow <- lapply(1:length(proteinMeansPerDilution), getCombIntensitiesMeans, 
                                       intensitiesPerDilution=proteinMeansPerDilution, intensityCombOrMean= "intensityProteinMean")
  proteinMeansPerDIAWorkflow.df <- do.call("rbind", proteinMeansPerDIAWorkflow)
  colnames(proteinMeansPerDIAWorkflow.df) <- c("proteinnames", "dia", "dilution", "intensity")
  proteinMeansPerDIAWorkflow.df <- addSoftwareLibraryColumns(proteinMeansPerDIAWorkflow.df, sepSpace = TRUE)
  proteinMeansPerDIAWorkflow.df
}

comInt.Ecoli <- getCombinedIntensitiesPerDIAWorkflow(diaWorkflowResults.log, humanOrEcoli="ECOLI")
comInt.Human <- getCombinedIntensitiesPerDIAWorkflow(diaWorkflowResults.log, humanOrEcoli="HUMAN")
comInt.EcoliAndHuman <- getCombinedIntensitiesPerDIAWorkflow(diaWorkflowResults.log, humanOrEcoli="")

protMeanInt.Ecoli <- getProteinMeansIntensitiesPerDIAWorkflow(diaWorkflowResults.log, humanOrEcoli="ECOLI")
protMeanInt.Human <- getProteinMeansIntensitiesPerDIAWorkflow(diaWorkflowResults.log, humanOrEcoli="HUMAN")
protMeanInt.EcoliAndHuman <- getProteinMeansIntensitiesPerDIAWorkflow(diaWorkflowResults.log, humanOrEcoli="")

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

give.AverageNumNonNAsPerSample <- function(x){
  return(c(y = max(x)*1.07, label = round(sum(!is.na(x))/23, digits = 2)))  # 23 is the number of samples per dilution
}

# SUPPLEMENTARY FIGURE S2
pdf(file="FigS2_comInt.Ecoli.pdf", width=10, height=8)
ggplot(comInt.Ecoli, aes(x = dilution, y = intensity)) + 
  geom_violin(scale="count") +
  geom_boxplot(width=0.08, outlier.size=0.5, outlier.alpha=0.4) +
  stat_summary(fun.data = give.AverageNumNonNAsPerSample, geom = "text", fun = max, size=2) +
  facet_grid(diaSoftware ~ diaLibrary, scales = "fixed") + 
  xlab("") +
  ylab("Log2(Intensity)") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(aes(yintercept = median(intensity, na.rm = TRUE)), colour = 'red') +
  ggtitle("E.coli")
dev.off()

pdf(file="diaWorkflowViolinPlot_comInt.Human.pdf", width=10, height=8)
ggplot(comInt.Human, aes(x = dilution, y = intensity)) + 
  geom_violin(scale="count") +
  geom_boxplot(width=0.08, outlier.size=0.5, outlier.alpha=0.4) +
  stat_summary(fun.data = give.AverageNumNonNAsPerSample, geom = "text", fun = max, size=2) +
  facet_grid(diaSoftware ~ diaLibrary, scales = "fixed") + 
  xlab("") +
  ylab("Log2(Intensity)") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(aes(yintercept = median(intensity, na.rm = TRUE)), colour = 'red') +
  ggtitle("Human")
dev.off()

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
for (i in 1:length(diaWorkflowResults.log)) {
  dia.name <- names(diaWorkflowResults.log)[i]
  df <- diaWorkflowResults.log[[i]]
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
combinedProteinNames <- lapply(diaWorkflowResults.log, function(x, combinedProteinNames){
  names <- unique(unlist(strsplit(row.names(x), ";")))
  combinedProteinNames <- c(combinedProteinNames, names)
  combinedProteinNames
}, combinedProteinNames=combinedProteinNames)

intersectProteinNames <- data.frame(table(unlist(combinedProteinNames)))
# Get protein names present in ALL DIA workflows
intersectProteinNamesAll <- intersectProteinNames[intersectProteinNames$Freq == 17,]$Var1

old <- c("CTGEF_HUMAN;CTGE9_HUMAN;CTGE8_HUMAN;CTGE6_HUMAN;CTGE4_HUMAN;MIA2_HUMAN",
         "IMA6_HUMAN;IMA7_HUMAN",
         "LIN7A_HUMAN;LIN7C_HUMAN",
         "NUD4B_HUMAN;NUDT4_HUMAN",
         "RRAGB_HUMAN;RRAGA_HUMAN",
         "SBNO1_HUMAN;SBNO2_HUMAN", 
         "WAC2A_HUMAN;WAC2C_HUMAN")
new <- c("MIA2_HUMAN",
         "IMA7_HUMAN",
         "LIN7C_HUMAN",
         "NUDT4_HUMAN",
         "RRAGA_HUMAN",
         "SBNO1_HUMAN",
         "WAC2A_HUMAN")


# blub <- list()
# for (i in 1:length(diaWorkflowResults.log)) {
#   df <- diaWorkflowResults.log[[i]]
#   df <- data.frame(protein.name = row.names(df))
# 
#   intersectBoolIntersect <- apply(df, 1, function(x) {
#     length(intersect(unlist(base::strsplit(x[1], ";")), intersectProteinNamesAll))>0
#   })
#  
#   for (j in 1:length(old)){
#     #print(df$protein.name == old[j])
#     try(df$protein.name[df$protein.name == old[j]] <- new[j])
#   }
#   
#   df.intersectProtNames <- data.frame(df[intersectBoolIntersect & df$protein.name != "LAP2A_HUMAN;LAP2B_HUMAN", ])
#   df.intersectProtNames <- df.intersectProtNames[order(df.intersectProtNames[,colnames(df.intersectProtNames)[1]]),]
#   blub <- list.append(blub, df.intersectProtNames)
#   
# }
# 
# blub.df <- do.call(cbind, blub)  
# colnames(blub.df) <- names(diaWorkflowResults.log)
# 
# uniques <- apply(blub.df, 1, function(x){
#   # sum(unique(x))
#   length(unique(x))
# })
# 
# blub.df <- data.frame(blub.df, uniques)
# final.intersect <- blub.df[,1]


diaWorkflowResults.log.intersect <- lapply(diaWorkflowResults.log, function(x){
  x <- x[row.names(x) !=  "LAP2A_HUMAN;LAP2B_HUMAN", ]
  for (j in 1:length(old)){
    #print(df$protein.name == old[j])
    try(row.names(x)[row.names(x) == old[j]] <- new[j])
  }
  
  x <- x[row.names(x) %in% intersectProteinNamesAll, ]
  x <- x[order(row.names(x)),]
  
  x
  # row.names(x)
})

# 
# diaWorkflowResults.log.intersect.df <- do.call(cbind, diaWorkflowResults.log.intersect)  
# colnames(diaWorkflowResults.log.intersect.df) <- names(diaWorkflowResults.log.intersect)
# 
# diaWorkflowResults.log.intersect.uniques <- apply(diaWorkflowResults.log.intersect.df, 1, function(x){
#   # sum(unique(x))
#   length(unique(x))
# })
# 
# diaWorkflowResults.log.intersect.df <- data.frame(diaWorkflowResults.log.intersect.df, diaWorkflowResults.log.intersect.uniques)


diaWorkflowResults.log.intersect <- lapply(diaWorkflowResults.log.intersect, function(x){
  unlist(x)
})

diaWorkflowResults.log.intersect.df <- do.call(cbind, diaWorkflowResults.log.intersect)  
colnames(diaWorkflowResults.log.intersect.df) <- names(diaWorkflowResults.log.intersect)


M <- cor(diaWorkflowResults.log.intersect.df, use="pairwise.complete.obs") # "pairwise.complete.obs" as there are NAs in regVals for some dias

pdf(file = paste0("corrplot_log2Intensities_diaworkflows.pdf"), width = 16, height = 16)
corrplot::corrplot.mixed(M,  upper = "ellipse", lower = "number",
                         tl.pos = "lt", tl.col = "black", mar=c(0,0,2,0))
dev.off()

pdf(file = paste0("networkplot_log2Intensities_diaworkflows.pdf"), width = 25, height = 25)
diaWorkflowResults.log.intersect.df %>% 
  corrr::correlate() %>% 
  corrr::network_plot()
dev.off()


########################################################################

plotUpsetPlot <- function(diaWorkflowResults.log, selectedDilutions=c("Lymph nodes", "Lymph nodes + 1:25 E.coli", 
                                                                      "Lymph nodes + 1:12 E.coli", "Lymph nodes + 1:6 E.coli"), str="") {
  diaWorkflowResults.log.names <- lapply(diaWorkflowResults.log, function(df){
    design.dilutions <-  addDilutionColumn(data.frame(dia=colnames(df)), colname="dia") 
    
    df <- df[, design.dilutions$dilution %in% selectedDilutions]
    # Remove empty rows
    df <- df[rowSums(is.na(df)) != ncol(df),]
    names <- unique(unlist(strsplit(row.names(df), ";")))
    names
  })
  
  # SUPPLEMENTARY FIGURE S3 AND S4
  
  lapply(c("E.coli & Human", "E.coli", "Human"), function(intersectionType, diaWorkflowResults.log.names) {
    print(intersectionType)
    if (intersectionType == "E.coli"){
      protein.names <- lapply(diaWorkflowResults.log.names, function(x){
        x[grepl("ECOLI", x)]
      })
    } else if (intersectionType == "Human"){
      protein.names <- lapply(diaWorkflowResults.log.names, function(x){
        x[grepl("HUMAN", x)]
      })
    } else if (intersectionType == "E.coli & Human"){
      protein.names <- diaWorkflowResults.log.names
    }
    
    pdf(paste0("FigS3AndS4_upsetPlot_", gsub(" ", "", intersectionType), "_noHEKIncluded_", str, ".pdf"), width=10, height=8)
    print(UpSetR::upset(fromList(protein.names),
                        order.by = "freq", 
                        nsets = length(protein.names),
                        sets = rev(names(protein.names)),
                        keep.order = TRUE, mainbar.y.label=paste0("Intersection Size ", intersectionType)))
    dev.off()
  }, diaWorkflowResults.log.names=diaWorkflowResults.log.names)
}


plotUpsetPlot(diaWorkflowResults.log, selectedDilutions=c("Lymph nodes", "Lymph nodes + 1:25 E.coli", 
                                                          "Lymph nodes + 1:12 E.coli", "Lymph nodes + 1:6 E.coli"), str="allDilutions") 

plotUpsetPlot(diaWorkflowResults.log, selectedDilutions=c("Lymph nodes + 1:25 E.coli", 
                                                          "Lymph nodes + 1:12 E.coli"), str="12and25Only") 



########################################################################

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

plotPCAGrid <- function(diaWorkflowResults.log, filename, labelled=FALSE) {
  pcaPlot.lst <- list()
  pcaPlot.lst <- lapply(1:length(diaWorkflowResults.log), function(idx){
    print(names(diaWorkflowResults.log)[idx])
    df <- diaWorkflowResults.log[[idx]]
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
      # ggtitle(names(diaWorkflowResults.log)[idx]) +
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
  
  p2_legend <- get_legend(pcaPlot.lst[[5]]) 
  
  grobArranged <- arrangeGrob(pcaPlot.lst[[5]] + theme(legend.position="none"), 
                              pcaPlot.lst[[2]] + theme(legend.position="none"), 
                              pcaPlot.lst[[3]] + theme(legend.position="none"),
                              pcaPlot.lst[[4]] + theme(legend.position="none"), 
                              pcaPlot.lst[[1]] + theme(legend.position="none"),
                              pcaPlot.lst[[12]] + theme(legend.position="none"), 
                              pcaPlot.lst[[9]] + theme(legend.position="none"), 
                              pcaPlot.lst[[10]] + theme(legend.position="none"),
                              pcaPlot.lst[[11]] + theme(legend.position="none"), 
                              patchwork::plot_spacer() + theme_minimal(), 
                              patchwork::plot_spacer() + theme_minimal(),
                              pcaPlot.lst[[6]] + theme(legend.position="none"), 
                              pcaPlot.lst[[7]] + theme(legend.position="none"), 
                              pcaPlot.lst[[8]] + theme(legend.position="none"), 
                              patchwork::plot_spacer() + theme_minimal(),
                              pcaPlot.lst[[17]] + theme(legend.position="none"), 
                              pcaPlot.lst[[13]] + theme(legend.position="none"), 
                              pcaPlot.lst[[15]] + theme(legend.position="none"),
                              pcaPlot.lst[[16]] + theme(legend.position="none"), 
                              pcaPlot.lst[[14]] + theme(legend.position="none"), 
                              nrow=4)
  
  
  g <- grid.arrange(rbind(tableGrob(t(c("PROSIT EDIA GPF", "DIANN AI GPF", "MaxQuant", "MSFragger", "Predicted")), theme = ttheme_minimal(), rows = ""), 
                          cbind(tableGrob(c("DIANN", "Skyline", "OSW", "Spectronaut"), theme = ttheme_minimal()), 
                                grobArranged,  size = "last"), size = "last"), p2_legend, 
                    nrow=2, heights=c(10, 2))
  
  ggsave(file=filename, g, width=15, height=12)
}
# dev.off()
#plotPCAGrid(diaWorkflowResults.log, filename="pcaNIPALS_allDIAWorkflows.pdf") 
#plotPCAGrid(diaWorkflowResults.log, filename="pcaNIPALS_allDIAWorkflows_labelled2.pdf", labelled=TRUE) 


diaWorkflowResults.log.normalized <- lapply(diaWorkflowResults.log, function(df){
  df <- df[, colnames(df) != "1-6_028"]
  df.norm <- limma::normalizeQuantiles(df)
  df.norm
})
#plotPCAGrid(diaWorkflowResults.log.normalized, filename="pcaNIPALS_allDIAWorkflows_QNnormalized.pdf") 

# SUPPLEMENTARY FIGURE S9
plotPCAGrid(diaWorkflowResults.log.normalized, filename="FigS9_pcaNIPALS_QNnormalized_labelled.pdf", labelled=TRUE) 



########################################################################

# GET NUMBER OF INTERSECT AND UNION PROTEINS FOR EACH BOOTSTRAP DATASET FOR EACH DIA WORKFLOW

#load("Converted_to_same_format_wide_plusIndices_plusIntersectingAndCombinedProteinNames.RData")
load("DIAsoftwareOutputProteinLevel_1to12And1to25Only_wideFormat_withBootstrapIndicesAndIntersectAndCombinedProteinNames.RData")

combinedProteinNames <- readRDS("combinedProteinNames.rds")
intersectProteinNames <- readRDS("intersectProteinNames.rds")

nEcoli.comb <- length(combinedProteinNames[grepl("_ECOLI", combinedProteinNames)])
# # [1] 2127
nHuman.comb <- length(combinedProteinNames[grepl("_HUMAN", combinedProteinNames)])
# # [1] 11516
# 
nEcoli.intersect <- length(intersectProteinNames[grepl("_ECOLI", intersectProteinNames)])
# # [1] 745
nHuman.intersect <- length(intersectProteinNames[grepl("_HUMAN", intersectProteinNames)])
# # [1] 4499


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


df[df$dia == "DIANN_DIANN_AI",]$dia <- "DIANN_Predicted"
df[df$dia == "Spectronaut_DirectDIA",]$dia <- "Spectronaut_Predicted"

df$dia <-gsub("_", " ", df$dia)
df <- df[order(df$dia),]
df <- df[,c(ncol(df),1:(ncol(df)-1))]

# Part of this table is depicted in SUPPLEMENTARY FIGURE S3 AND S4
write.csv(df, file = "numberOfProteinsIntersectCombined.csv", row.names = FALSE)


#####################################################################################################################
#####################################################################################################################
# BOOTSTRAP

# READ IN FILES AND FORMAT RESULT DATA FRAME
theme_set(theme_minimal())

listfiles <- list.files(path = "csvs", full.names = TRUE, pattern = "*.csv") 

listfiles.entrynumbers <- lapply(X = listfiles, FUN = function(x) {
  length(count.fields(x, skip = 1))
})
names(listfiles.entrynumbers) <- listfiles
listfiles.entrynumbers <- t(data.frame(listfiles.entrynumbers))

df.combined <- ldply(listfiles, readr::read_csv)
colnames(df.combined)[which(names(df.combined) == "sparcityReduction")] <- "sparsityReduction"
colnames(df.combined)[which(names(df.combined) == "sparcRed.runtime")] <- "sparsRed.runtime"

df.combined[df.combined$X1 == 2,]$X1 <- 1

df.combined[(df.combined$dia == "DIANN_DIANN_AI_GPF") 
            & (df.combined$normalization == "median") & (df.combined$sparsityReduction == "SR66") 
            & (df.combined$statTest == "LM") & (df.combined$groupSize == 4) & (df.combined$X1 == 25201), ]$X1 <- 2

# The entry of row 21000 in first columns of _3 and _6.csvs has been changed manually to 210000
df.combined[df.combined$X1 == 210000,]$X1 <- 2
df.combined$X1 <- df.combined$X1 %% 2100
df.combined[df.combined$X1 == 0,]$X1 <- 2100

factor.cols <- c("dia", "normalization", "sparsityReduction", "statTest")
df.combined[,factor.cols] <- lapply(df.combined[,factor.cols], factor) 

lapply(df.combined, class)

########################################################################
# Add additional bootstrap data characteristics 

dia.names <- as.character(sort(unique(df.combined$dia)))

datalist <- list()

for (dia.name in dia.names) {
  df <- read.csv(paste0("20210625_additionalCharacteristics/", dia.name, "_additionalCharactInfos_allBootstrapDatasets.csv"))
  df$dia <- dia.name
  names(df)[names(df) == 'X'] <- 'X1'
  datalist <- list.append(datalist, df)
}

df.addInfo <- do.call(rbind, datalist)

# Match dia and X1
df.combined <- df.combined %>% inner_join(df.addInfo, by=c("dia","X1"))

df.combined$percNATotal2 <- NULL

##########################################################################
# Add pAUCs
pAUC.csvs <- list.files(path = "pAUC_statsdfAdditionalInfos", full.names = TRUE, pattern = "*.csv") 

listfiles.entrynumbers <- lapply(X = pAUC.csvs, FUN = function(x) {
  length(count.fields(x, skip = 1))
})
names(listfiles.entrynumbers) <- pAUC.csvs
listfiles.entrynumbers <- t(data.frame(listfiles.entrynumbers))

df.combined.pAUCs <- plyr::ldply(pAUC.csvs, readr::read_csv)
df.combined <- df.combined %>% inner_join(df.combined.pAUCs, by=c("dia", "X1", "normalization", "sparsityReduction", 
                                                                  "statTest",  "groupSize"))

##########################################################################
RMSE.csvs <- list.files(path = "correctedRMSE", full.names = TRUE, pattern = "*.csv") 

listfiles.entrynumbers <- lapply(X = RMSE.csvs, FUN = function(x) {
  length(count.fields(x, skip = 1))
})
names(listfiles.entrynumbers) <- RMSE.csvs
listfiles.entrynumbers <- t(data.frame(listfiles.entrynumbers))

df.combined.RMSE <- plyr::ldply(RMSE.csvs, readr::read_csv)
df.combined <- df.combined %>% inner_join(df.combined.RMSE, by=c("dia", "X1", "normalization", "sparsityReduction", 
                                                                 "statTest",  "groupSize"))
##########################################################################
`%notin%` <- Negate(`%in%`)

# Remove LM and Lasso
df.combined <- df.combined[df.combined$statTest  %notin% c("LM", "Lasso"),]

df.combined$normalization <- factor(df.combined$normalization, levels = c("unnormalized","median", "QN", "TRQN"))
df.combined$sparsityReduction <- factor(df.combined$sparsityReduction, levels = c("NoSR", "SR66", "SR90"))
df.combined$statTest <- factor(df.combined$statTest, levels = c("ttest", "GLMgamma", "limma", "ROTS", "SAM")) # Lasso and LM removed
df.combined$log1pSparsRed.runtime <- log1p(df.combined$sparsRed.runtime)
df.combined$log1pNormalization.runtime <- log1p(df.combined$normalization.runtime)
df.combined$log1pStatTest.runtime <- log1p(df.combined$statTest.runtime)
df.combined$portionEcoliProteins <- df.combined$nEcoliProteins/df.combined$nAllProteins
df.combined$nAllProteins.pre <- df.combined$nEcoliProteins.pre + df.combined$nHumanProteins.pre
df.combined$portionEcoliProteins.pre <- df.combined$nEcoliProteins.pre/df.combined$nAllProteins.pre


# Remove the following measures as they are deprecated 
df.combined <- df.combined[, colnames(df.combined) %notin% c("p.auc", "p.auk", "p.PrecRecallauc", "p.auc.pre", 
                                                             "p.auk.pre", "p.PrecRecallauc.pre", "p.auc.comb", "p.auk.comb", 
                                                             "p.PrecRecallauc.comb", "p.auc.intersect", "p.auk.intersect", 
                                                             "p.PrecRecallauc.intersect", "fc.auc", "fc.auc.intersectProtNames", 
                                                             "fc.auc.FCpValRanksum", "fc.auc.FCpValRanksum.intersectProtNames",
                                                             "RMSEEcoli", "RMSEHuman", "RMSEHumanAndEcoli",
                                                             "RMSEEcoli.Intersect", "RMSEHuman.Intersect",
                                                             "RMSEHumanAndEcoli.Intersect",
                                                             "NRMSEEcoli", "NRMSEHuman", "NRMSEHumanAndEcoli",
                                                             "NRMSEEcoli.Intersect", "NRMSEHuman.Intersect",
                                                             "NRMSEHumanAndEcoli.Intersect",
                                                             "p.pauc_0_correctFALSE.combined",
                                                             "p.pauc_0_correctFALSE.intersect",
                                                             "p.pauc_0_correctFALSE.DiaWorkflowProteins")] 


df.combined[df.combined$dia == "DIANN_DIANN_AI",]$dia <- "DIANN_Predicted"
df.combined[df.combined$dia == "Spectronaut_DirectDIA",]$dia <- "Spectronaut_Predicted"
df.combined$dia <- gsub("_", " ", df.combined$dia)

#write.csv(df.combined, "20210704_benchmarkCombined.csv", row.names = FALSE)
#saveRDS(df.combined, "20210704_benchmarkCombined.rds")

settings <- c("dia", "normalization", "sparsityReduction", "statTest") 

data.characteristics <- c("groupSize", "portionEcoliProteins", "nAllProteins", 
                          "portionEcoliProteins.pre", "nAllProteins.pre", 
                          "medianSampleVariance", "medianProteinVariance", "KS.SignProp", 
                          "percNATotal", "percOfRowsWithNAs",
                          "entropy.wNAs", "kurtosis.wNAs", "meanDeviation.wNAs", 
                          "skewness.wNAs", "uniformity.wNAs", "variance.wNAs", "RMS.wNAs", 
                          "var.groups.ratio.wNAs", "entropy.woNAs", "kurtosis.woNAs", "meanDeviation.woNAs", 
                          "skewness.woNAs", "uniformity.woNAs", "variance.woNAs", "RMS.woNAs", 
                          "var.groups.ratio.woNAs", "prctPC1.woNAs", "elongation.woNAs", 
                          "flatness.woNAs", "heterosc.oneMinuspi0", "nProteins.woNAs")

eval.measures <- c("p.pauc_0.8_correctFALSE.combined", 
                   "p.pauc_0.9_correctFALSE.combined", "p.pauc_0.95_correctFALSE.combined", 
                   "sensAtpVal005.combined", "regpValsProp.combined", 
                   "rankAUC.FCpVal.combined", "rankAUC.FC.combined", "rankAUC.pVal.combined", 
                   "p.pauc_0.8_correctFALSE.intersect", 
                   "p.pauc_0.9_correctFALSE.intersect", "p.pauc_0.95_correctFALSE.intersect", 
                   "sensAtpVal005.intersect", 
                   "regpValsProp.intersect", "rankAUC.FCpVal.intersect", "rankAUC.FC.intersect", 
                   "rankAUC.pVal.intersect", 
                   "p.pauc_0.8_correctFALSE.DiaWorkflowProteins", "p.pauc_0.9_correctFALSE.DiaWorkflowProteins", 
                   "p.pauc_0.95_correctFALSE.DiaWorkflowProteins", 
                   "sensAtpVal005.DiaWorkflowProteins", "regpValsProp.DiaWorkflowProteins", 
                   "rankAUC.FCpVal.DiaWorkflowProteins", "rankAUC.FC.DiaWorkflowProteins", 
                   "rankAUC.pVal.DiaWorkflowProteins",
                   "RMSEEcoli.DiaWorkflowProteins", "RMSEHuman.DiaWorkflowProteins", "RMSEHumanAndEcoli.DiaWorkflowProteins",
                   "RMSEEcoli.intersect", "RMSEHuman.intersect",
                   "RMSEHumanAndEcoli.intersect"
)

eval.measure <- "p.pauc_0.8_correctFALSE.DiaWorkflowProteins"

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


ggsave(filename="violinPlots_settingsCombined_statTest_statTest.runtime_horizontal.pdf", plot=plotViolinPlotsCombSettingHorizontal("statTest", df.combined, "statTest.runtime"), height = 3, width = 5)

# FIGURE 4

eval.measures <- c("p.pauc_0.8_correctFALSE.DiaWorkflowProteins", "p.pauc_0.8_correctFALSE.intersect", "p.pauc_0.8_correctFALSE.combined")

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

eval.measure <- "p.pauc_0.8_correctFALSE.DiaWorkflowProteins"
selectedSR <- "NoSR"

df.combined.SelectedSRonly <- df.combined[df.combined$sparsityReduction==selectedSR,]
plotsSettings.SelectedSR <- lapply(setdiff(settings, "sparsityReduction"), plotViolinPlotsCombSetting,df.combined=df.combined.SelectedSRonly, eval.measure=eval.measure)

pdf(file = paste0("violinPlots_settingsCombined_", selectedSR, "_", eval.measure,".pdf"), height = 18, width = 5)
cowplot::plot_grid(plotlist=plotsSettings.SelectedSR, align = "v", ncol=1)
dev.off()

eval.measure <- "p.pauc_0.8_correctFALSE.DiaWorkflowProteins"
pdf(file = paste0("Fig4B_violinPlots_settingsCombined_", eval.measure, "_SRs.pdf"), height = 5, width = 5)
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
   ggplot2::ggsave(file = paste0(str, eval.measure, "_", settChar,".pdf"), plot = ggplot.violin, device = "pdf", height = 8, width = 14)
}

df.combined.LibSoftware <- df.combined
df.combined.LibSoftware <- addSoftwareLibraryColumns(df.combined.LibSoftware, sepSpace=TRUE)


# lapply(eval.measures, function(eval.measure){
#   lapply(settChars, getViolinGrid, eval.measure = eval.measure, df=df.combined.LibSoftware, str="dfcombinedsel_facet_violin_")
# })

df.combined.LibSoftware.NoSROnly <- df.combined.LibSoftware[df.combined.LibSoftware$sparsityReduction=="NoSR",]
# lapply(eval.measures, function(eval.measure){
#   lapply(settChars, getViolinGrid, eval.measure = eval.measure, df=df.combined.LibSoftware.NoSROnly, str="dfcombinedsel_facet_violin_NoSROnly_")
# })

# SUPPLEMENTARY FIGURE S6
getViolinGrid("normalization","p.pauc_0.8_correctFALSE.DiaWorkflowProteins", df.combined.LibSoftware.NoSROnly, str="FigS6_normalization_violin_NoSROnly_")

# SUPPLEMENTARY FIGURE S7
getViolinGrid("statTest","p.pauc_0.8_correctFALSE.DiaWorkflowProteins", df.combined.LibSoftware.NoSROnly, str="FigS7_statTest_violin_NoSROnly_")


# getViolinGridSplitBySR <- function(settChar, eval.measure, df, str){
#   ggplot.violin <- ggplot(df, aes(x = get(settChar), y = get(eval.measure), fill=sparsityReduction)) + 
#     xlab(settChar) + ylab(eval.measure) +
#     geom_boxplot(outlier.size=0.1, outlier.alpha=0.4, alpha=0.5) +
#     ggplot2::scale_fill_brewer(palette="Blues")  +
#     facet_grid(diaSoftware ~ diaLibrary, scales = "fixed") +
#     theme_bw() +
#     xlab("") +
#     theme(strip.placement = "outside", axis.text.x = element_text(angle = 55, hjust = 1)) +
#     geom_hline(aes(yintercept = median(get(eval.measure), na.rm = TRUE)), colour = 'red') 
#   ggplot2::ggsave(file = paste0(str, eval.measure, "_", settChar,".pdf"), plot = ggplot.violin, device = "pdf", height = 8, width = 14)
# }
# 
# evals <-c("p.pauc_0.8_correctFALSE.DiaWorkflowProteins", "p.pauc_0.8_correctFALSE.intersect", "p.pauc_0.8_correctFALSE.combined")
# for (eval in evals){
#   getViolinGridSplitBySR(settChar="statTest", eval.measure = eval, df=df.combined.LibSoftware, str="dfcombinedsel_facet_violin_separatedBySR_")
# }


df.combined.LibSoftware.NoSROnly$groupSize <- as.numeric(as.character(df.combined.LibSoftware.NoSROnly$groupSize))

# SUPPLEMENTARY FIGURE S11
for (eval.measure2 in c("p.pauc_0.8_correctFALSE.DiaWorkflowProteins", 
                        "p.pauc_0.8_correctFALSE.intersect", 
                        "p.pauc_0.8_correctFALSE.combined",
                        "sensAtpVal005.DiaWorkflowProteins", 
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
  ggsave(paste0("FigS11_ggplot_medianLine_groupSize_NoSROnly_", eval.measure2, ".pdf"), ggplot.medianLine.groupSize.statTest, height = 8, width = 10)
}

df.combined.LibSoftware.selected <- df.combined.LibSoftware[, c("dia",
                                                                "sparsityReduction",
                                                                "groupSize",
                                                                "statTest",
                                                               "p.pauc_0.8_correctFALSE.combined",
                                                               "p.pauc_0.8_correctFALSE.intersect",
                                                                "p.pauc_0.8_correctFALSE.DiaWorkflowProteins")]

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

  ggsave(paste0("medianLine_SRProteinSetstatTest_", gsub(" ", "_", selectedDia), ".pdf"), ggplot.medianLine.groupSize.SRProteinSetstatTest, height = 8, width = 10)
}

########################################################################

data.characteristics.sel <- c("groupSize",  
                              "medianSampleVariance", "medianProteinVariance", 
                              "percNATotal", "kurtosis.wNAs", "skewness.wNAs", 
                              "var.groups.ratio.wNAs", "prctPC1.woNAs")

eval.measures.sel <- sort(c("p.pauc_0.8_correctFALSE.combined", 
                            "sensAtpVal005.combined", "regpValsProp.combined", 
                            "p.pauc_0.8_correctFALSE.intersect",
                            "sensAtpVal005.intersect", "regpValsProp.intersect", 
                            "p.pauc_0.8_correctFALSE.DiaWorkflowProteins", 
                            "sensAtpVal005.DiaWorkflowProteins", "regpValsProp.DiaWorkflowProteins", 
                            "RMSEEcoli.DiaWorkflowProteins", "RMSEHuman.DiaWorkflowProteins", 
                            "RMSEEcoli.intersect", "RMSEHuman.intersect"))

sett <- "dia"  
sel <- c(data.characteristics.sel, eval.measures.sel)

df.combined.pairsplotsel <- df.combined[df.combined$sparsityReduction=="NoSR", c(sett, sel)]
df.combined.pairsplotsel$groupSize <- as.numeric(as.character(df.combined.pairsplotsel$groupSize))

# SUPPLEMENTARY FIGURE S10
for (dia in sort(unique(df.combined.pairsplotsel$dia))){
  df.combined.pairsplotsel.dia <- df.combined.pairsplotsel[df.combined.pairsplotsel$dia==dia,]
  df.combined.pairsplotsel.dia$dia <- NULL
  
  M <- cor(df.combined.pairsplotsel.dia, use="pairwise.complete.obs") # "pairwise.complete.obs" as there are NAs in regVals for some dias
  
  pdf(file = paste0("FigS10_corrplot_charact_NoSR_", dia, ".pdf"), width = 15, height = 15)
  corrplot.mixed(M,  upper = "ellipse", lower = "number",
                 tl.pos = "lt", tl.col = "black", main=dia, mar=c(0,0,2,0))
  dev.off()
  
  # pdf(file = paste0("corrplot2_charact_NoSR_", dia, ".pdf"), width = 15, height = 15)
  # corrplot(M, method = "number", tl.col="black", main=dia, mar=c(0,0,2,0))
  # dev.off()
}

df.combined.pairsplotsel.long <- reshape2::melt(df.combined.pairsplotsel[, c("dia", setdiff(data.characteristics.sel, "groupSize"))])
df.combined.pairsplotsel.long <- addSoftwareLibraryColumns(df.combined.pairsplotsel.long, sepSpace = TRUE)

# SUPPLEMENTARY FIGURE S5
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
ggsave(file=paste0("FigS5_charact_bootstrapDatasets_NoSROnly.pdf"), ggplot.charact, height=14, width=11)

########################################################################

sett <- "dia"
sel <- c("X1", "groupSize", settings)

eval.measures <- c("p.pauc_0.8_correctFALSE.DiaWorkflowProteins", "p.pauc_0.8_correctFALSE.intersect", "p.pauc_0.8_correctFALSE.combined")
for (eval.measure in eval.measures){
  df.combined.eval <- df.combined[, c(sel, eval.measure)]
  
  df.combined.eval$concatenated <- do.call(paste, c(df.combined.eval[, setdiff(sel, sett)], sep = "_"))
  df.combined.eval <- df.combined.eval[, c(sett, "concatenated", eval.measure)]
  
  df.combined.eval.wide <- spread(df.combined.eval, key = concatenated, value = get(eval.measure))
  row.names(df.combined.eval.wide) <- df.combined.eval.wide[,sett]
  df.combined.eval.wide[,sett] <- NULL
  df.combined.eval.wide <- t(df.combined.eval.wide)
  
  M <- cor(df.combined.eval.wide, use="pairwise.complete.obs")
  
  # SUPPLEMENTARY FIGURE S8
  pdf(file = paste0("FigS8_corrplot_", sett, "_",eval.measure,".pdf"), width = 15, height = 15)
  corrplot.mixed(M,  upper = "ellipse", lower = "number",
                 tl.pos = "lt", tl.col = "black")
  dev.off()
  
  pdf(file = paste0("corrplot2_", sett, "_",eval.measure,".pdf"), width = 15, height = 15)
  corrplot(M, method = "number", tl.col="black")
  dev.off()
  
  pdf(file = paste0("networkplot_", sett, "_",eval.measure,".pdf"), width = 20, height = 20)
  print(df.combined.eval.wide %>% 
    corrr::correlate() %>% 
    corrr::network_plot())
  dev.off()
}

session <- sessionInfo()
sink("sessionInfo_visualizations.txt")
print(session)
sink()

########################################################################
# Ranking

library(dplyr)


df.combined2 <- df.combined
df.combined2 <- addSoftwareLibraryColumns(df.combined2, sepSpace = TRUE)

diaSpecifications <- c("dia", "diaSoftware", "diaLibrary")
settings <- c("statTest", "sparsityReduction", "normalization")
proteinLists <- c("DiaWorkflowProteins", "intersect", "combined")
# ranks.statTest <- df.combined2 %>% 
#   dplyr::group_by(X1, statTest, dia) %>% 
#   summarise(mean.p.pauc_0.8_correctFALSE.DiaWorkflowProteins = mean(p.pauc_0.8_correctFALSE.DiaWorkflowProteins)) %>% 
#   dplyr::group_by(X1, dia) %>%
#   mutate(ranking = rank(- mean.p.pauc_0.8_correctFALSE.DiaWorkflowProteins)) %>% 
#   dplyr::group_by(statTest, dia) %>% 
#   summarise(mean.ranking = mean(ranking))


for (setting in settings){
  for (diaSpecification in diaSpecifications){
    bump.lst <- list()
    for (proteinList in proteinLists) {
      ranks.df <- df.combined2 %>% 
        dplyr::group_by(X1, !!as.name(setting), dia, diaSoftware, diaLibrary) %>% 
        summarise(mean.pAUC = mean(!!as.name(paste0("p.pauc_0.8_correctFALSE.", proteinList)))) %>% 
        dplyr::group_by(X1, dia) %>%
        mutate(ranking = rank(- mean.pAUC)) %>% 
        dplyr::group_by(!!as.name(setting), !!as.name(diaSpecification)) %>% 
        summarise(mean.ranking = mean(ranking))
      
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
        theme(legend.title = element_blank()) +
        # xlab(setting) +
        xlab("") +
        scale_alpha(guide = 'none') +
        guides(colour = guide_legend(override.aes = list(alpha = 0.5)))
      ggsave(paste0("bumpchart_", diaSpecification, "_", setting, "_", proteinList, ".pdf"), bump.chart, width=width, height=5)
      
    }
    
    meanOverProteinLists.df <- data.frame(bump.lst[[1]][,1:2], mean.mean.ranking=rowMeans(data.frame(bump.lst[[1]][,"mean.ranking"], bump.lst[[2]][,"mean.ranking"], bump.lst[[3]][,"mean.ranking"])))
    bump.chart <- ggplot(data = meanOverProteinLists.df, aes(x = get(setting), y = mean.mean.ranking, group = get(diaSpecification), alpha=0.5)) +
      geom_line(aes(color = get(diaSpecification)), size = 1) +
      geom_point(aes(color = get(diaSpecification)), size = 2) +
      scale_y_reverse(limits=c(length(unique(ranks.df[[setting]],0)), 1), breaks = 1:nrow(meanOverProteinLists.df)) +
      ggtitle(paste0(diaSpecification, ", ", setting)) +
      theme(legend.title = element_blank()) +
      # xlab(setting) +
      xlab("") +
      scale_alpha(guide = 'none') +
      guides(colour = guide_legend(override.aes = list(alpha = 0.5)))
    ggsave(paste0("bumpchart_", diaSpecification, "_", setting, "_meanOverproteinLists.pdf"), bump.chart, width=width, height=5)
  }
}

