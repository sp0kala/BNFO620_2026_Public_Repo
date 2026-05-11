install.packages("remotes")
remotes::install_github("Mordziarz/RSCUcaller")

library(RSCUcaller)

install.packages(setdiff(c("stats", "dplyr", "ggplot2", 
                           "ggpubr", "seqinr","rstatix",
                           "patchwork", "forcats", "phylogram",
                           "circlize", "smplot2", "stringr"), 
                         installed.packages()[,"Package"]))

lapply(c("stats", "dplyr", "ggplot2", 
         "ggpubr", "seqinr","rstatix",
         "patchwork", "forcats", "phylogram",
         "circlize", "smplot2", "stringr"), library, character.only = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"):
  
  BiocManager::install(setdiff(c("ComplexHeatmap","ggtree"), 
                               installed.packages()[,"Package"]))

lapply(c("ComplexHeatmap","ggtree"), 
       library, character.only = TRUE)

library(ggplot2)
library(dplyr)
library(seqinr)
library(rstatix)
library(ggpubr)
library(patchwork)
library(forcats)
library(phylogram)
library(circlize)
library(ComplexHeatmap)
library(smplot2)
library(ggtree)
library(stats)
library(stringr)

install.packages("devtools")
library(devtools)
devtools::install_github('Mordziarz/RSCUcaller')
library(RSCUcaller)

set.seed(123)

getwd()
setwd("C:/Users/User/OneDrive/Desktop/mitogenomes")

getwd()

path1 <- "T_acuta_annotated.fasta"
path2 <- "T_alcocki_annotated.fasta"
path3 <- "T_arcuata_annotated.fasta"
path4 <- "T_capricornis_annotated.fasta"
path5 <- "T_paradussumieri_annotated.fasta"
path6 <- "T_polita_annotated.fasta"
path7 <- "T_rosea_annotated.fasta"

samples_table <- data.frame(
  sequence_path = c(path1, path2, path3, path4, path5, path6, path7),
  sample_name = c("1_T_acuta", "2_T_alcocki", "3_T_arcuata", "4_T_capricornis", 
                  "5_T_paradussumieri", "6_T_polita", "7_T_rosea")
)

prepare_fasta(
  samples_table = samples_table,
  file_out = "your_fasta.fasta"
)

get_RSCU_out <- get_RSCU(merged_sequences = "your_fasta.fasta")

get_RSCU_other(merged_sequences = "your_fasta.fasta",codon_table_id = 5,pseudo_count = 0)

get_matrix(get_RSCU_out = get_RSCU_out)

heatmap_RSCU(get_RSCU_out = get_RSCU_out, select = "heatmap", heatmap_color = "red_blue")

single_histogram(get_RSCU_out = get_RSCU_out, title = expression(italic("#Insert Name Here#")))

