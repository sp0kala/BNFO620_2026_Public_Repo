#
# Program to concatenate aligned gene fasta file
#

library(seqinr)

setwd("/Users/stellacastro/Downloads/BNFO620/BNFO620_2026_3/Aligned_Sequences/genes/aligned_genes")
files <- c("atp6_aligned.fas", "atp8_aligned.fas", "cob_aligned.fas", "cox1_aligned.fas",
           "cox2_aligned.fas", "cox3_aligned.fas", "nad1_aligned.fas", "nad2_aligned.fas", "nad3_aligned.fas",
           "nad4_aligned.fas", "nad4l_aligned.fas", "nad5_aligned.fas", "nad6_aligned.fas")

# read in all files with seqinr
genes <- lapply(files, read.fasta, seqtype = "DNA", as.string = TRUE)
# make sure all have the same header (ex: >T_Capricornis_cox1 to just >T_Capriconris)
label_headers <- function(gene){
  names(gene) <- sub("_[^_]+$", "", names(gene))
  return(gene)
}
genes <- lapply(genes, label_headers)

# save species names
species <- names(genes[[1]])


# list to store new fasta info
concatenated <- list()

# go through each species
for (s in species){
  combined_seq <- ""
  # go through each gene
  for (g in genes){
    # concatenate new gene wiith combined_seq
    combined_seq <- paste0(combined_seq, g[[s]])
  }
  # add species to list
  concatenated[[s]] <- combined_seq
}

# write to fasta file
write.fasta(
  sequences = concatenated,
  names = names(concatenated),
  file.out = "concatenated.fasta",
  as.string = TRUE
)
