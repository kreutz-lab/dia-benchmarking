
#load libraries
library(tidyverse)


#set working directory to current folder of R script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Skyline input
Skyline_MaxQuant      <- readRDS("MaxQuant_msstats_noNorm_noImpu.rda")
Skyline_MSFragger     <- readRDS("msfragger_noNorm_noImpu.rda")
Skyline_DIANN_AI_GPF  <- readRDS("DIANN_GPF_msstats_noNorm_noImpu.rda")
Skyline_PROSIT_EDIA_GPF <- readRDS("PROSIT_GPF_noNorm_noImpu.rda")

# run level data "Run", "Protein.Names", "Species" , "PG.Quantity"
RunLevel_MaxQuant <- Skyline_MaxQuant %>% 
  select(originalRUN, Protein, LogIntensities) %>%
  mutate(Protein.Names = as.character(Protein), Protein = NULL) %>%
  rename(Run = originalRUN, PG.Quantity = LogIntensities) %>%
  mutate(Species = substr(.$Protein.Names, start = nchar(Protein.Names)-4, stop = nchar(Protein.Names))) %>%
  mutate_if(is.numeric, function(x) 2^x)


RunLevel_MSFragger <- Skyline_MSFragger %>% 
  select(originalRUN, Protein, LogIntensities) %>%
  mutate(Protein.Names = as.character(Protein), Protein = NULL) %>%
  rename(Run = originalRUN, PG.Quantity = LogIntensities) %>%
  mutate(Species = substr(.$Protein.Names, start = nchar(Protein.Names)-4, stop = nchar(Protein.Names))) %>%
  mutate_if(is.numeric, function(x) 2^x)


RunLevel_DIANN_AI_GPF <- Skyline_DIANN_AI_GPF %>% 
  select(originalRUN, Protein, LogIntensities) %>%
  mutate(Protein.Names = as.character(Protein), Protein = NULL) %>%
  rename(Run = originalRUN, PG.Quantity = LogIntensities) %>%
  mutate(Species = substr(.$Protein.Names, start = nchar(Protein.Names)-4, stop = nchar(Protein.Names))) %>%
  mutate_if(is.numeric, function(x) 2^x)


RunLevel_PROSIT_EDIA_GPF <- Skyline_PROSIT_EDIA_GPF %>% 
  select(originalRUN, Protein, LogIntensities) %>%
  mutate(Protein.Names = as.character(Protein), Protein = NULL) %>%
  rename(Run = originalRUN, PG.Quantity = LogIntensities) %>%
  mutate(Species = substr(.$Protein.Names, start = nchar(Protein.Names)-4, stop = nchar(Protein.Names))) %>%
  mutate_if(is.numeric, function(x) 2^x)



rm(Skyline_DIANN_AI_GPF)
rm(Skyline_MaxQuant)
rm(Skyline_MSFragger)
rm(Skyline_PROSIT_EDIA_GPF)



save.image("Converted_to_same_format_Skyline_NoNorm_NoImp.RData")

write.csv(RunLevel_DIANN_AI_GPF, "RunLevelDIANN.csv")



