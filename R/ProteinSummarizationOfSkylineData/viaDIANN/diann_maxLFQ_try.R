#load libraries
library(tidyverse)
library(diann)


#set working directory to current folder of R script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

DIANN_results <- read.csv("SkylineNEWresults/Skyline_DIANN_GPF_msstats.csv", na.strings = "#N/A") 

DIANN_processed <- DIANN_results %>%
  filter(Detection.Q.Value < 0.01) %>%
  filter(Truncated == "False") %>%
  mutate(Precursor.Id = paste(Peptide.Modified.Sequence, Precursor.Charge, Fragment.Ion, sep = "_"))
  


proteins.maxlfq <- diann_maxlfq(x = DIANN_processed, 
                                sample.header = "File.Name",
                                group.header="Protein.Name", 
                                id.header = "Precursor.Id", 
                                quantity.header = "Area")


PROSIT.maxlfq <- read.csv("SkylineNEWresults/Skyline_PROSIT_GPF_msstats.csv", na.strings = "#N/A") %>%
  filter(Detection.Q.Value < 0.01) %>%
  filter(Truncated == "False") %>%
  mutate(Precursor.Id = paste(Peptide.Modified.Sequence, Precursor.Charge, Fragment.Ion, sep = "_")) %>%
  diann_maxlfq(           sample.header = "File.Name",
                                group.header="Protein.Name", 
                                id.header = "Precursor.Id", 
                                quantity.header = "Area")


MSFragger.maxlfq <- read.csv("SkylineNEWresults/Skyline_msfragger_msstats.csv", na.strings = "#N/A") %>%
  filter(Detection.Q.Value < 0.01) %>%
  filter(Truncated == "False") %>%
  mutate(Precursor.Id = paste(Peptide.Modified.Sequence, Precursor.Charge, Fragment.Ion, sep = "_")) %>%
  diann_maxlfq(           sample.header = "File.Name",
                          group.header="Protein.Name", 
                          id.header = "Precursor.Id", 
                          quantity.header = "Area")

MaxQuant.maxlfq <- read.csv("SkylineNEWresults/Skyline_MaxQuant_msstats.csv", na.strings = "#N/A") %>%
  filter(Detection.Q.Value < 0.01) %>%
  filter(Truncated == "False") %>%
  mutate(Precursor.Id = paste(Peptide.Modified.Sequence, Precursor.Charge, Fragment.Ion, sep = "_")) %>%
  diann_maxlfq(           sample.header = "File.Name",
                          group.header="Protein.Name", 
                          id.header = "Precursor.Id", 
                          quantity.header = "Area")


DIANN_export <- proteins.maxlfq %>% 
  as.data.frame() %>%
  mutate(Protein.Names = rownames(.)) %>%
  pivot_longer(cols = -Protein.Names, names_to = "Run", values_to = "PG.Quantity") %>%
  mutate(Species = substr(.$Protein.Names, start = nchar(Protein.Names) - 4, stop = nchar(Protein.Names)))

PROSIT_GPF_export <- PROSIT.maxlfq %>% 
  as.data.frame() %>%
  mutate(Protein.Names = rownames(.)) %>%
  pivot_longer(cols = -Protein.Names, names_to = "Run", values_to = "PG.Quantity") %>%
  mutate(Species = substr(.$Protein.Names, start = nchar(Protein.Names) - 4, stop = nchar(Protein.Names)))

MSFragger_export <- MSFragger.maxlfq %>% 
  as.data.frame() %>%
  mutate(Protein.Names = rownames(.)) %>%
  pivot_longer(cols = -Protein.Names, names_to = "Run", values_to = "PG.Quantity") %>%
  mutate(Species = substr(.$Protein.Names, start = nchar(Protein.Names) - 4, stop = nchar(Protein.Names)))

MaxQuant_export <- MaxQuant.maxlfq %>% 
  as.data.frame() %>%
  mutate(Protein.Names = rownames(.)) %>%
  pivot_longer(cols = -Protein.Names, names_to = "Run", values_to = "PG.Quantity") %>%
  mutate(Species = substr(.$Protein.Names, start = nchar(Protein.Names) - 4, stop = nchar(Protein.Names)))





