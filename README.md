Comparative Analysis of Mitochondrial Genomes of the Tubuca Genus

A comparative bioinformatics study analyzing mitochondrial genome structure, codon usage bias, natural selection, and phylogenetic relationships across seven Tubuca species (fiddler crabs).

Overview

This project investigates the mitochondrial genomes (mitogenomes) of seven Tubuca species using comparative genomics and phylogenetic analysis techniques. The study combines genome annotation, codon usage analysis, Ka/Ks substitution analysis, and phylogenetic tree construction to better understand evolutionary relationships within the genus.

Species Included
Tubuca acuta
Tubuca arcuata
Tubuca capricornis
Tubuca paradussumieri
Tubuca polita
Tubuca rosea
Tubuca alcocki

Outgroup:
Gelasimus borealis

Project Goals
Characterize mitochondrial genome organization and composition
Analyze codon usage bias using RSCU values
Evaluate synonymous vs. nonsynonymous substitution rates (Ka/Ks)
Construct phylogenetic relationships using:
Maximum Likelihood (MEGA)
Bayesian inference (RevBayes)
Methods
Genome Assembly & Annotation
Whole Genome Sequencing performed on Illumina HiSeqX

Annotation performed using:
MitoZ
MITOS2
Galaxy Server
Comparative Genomics
Genome visualization with Proksee
GC/AT skew calculations
Intergenic region analysis using BLAST
Codon Usage Analysis
RSCU values generated with RSCUcaller
Heatmaps and dendrograms generated in R using ggplot2
Selection Analysis

Ka/Ks ratios computed in MEGA using:
Nei-Gojobori method
Jukes-Cantor correction
Phylogenetic Analysis

Gene extraction and alignment:
MUSCLE alignment
Codon-based alignment
Tree construction:
Maximum Likelihood (GTR model)
Bayesian MCMC analysis in RevBayes
Key Findings
Genome Conservation

Genome sizes ranged from:
15,614 bp – 15,955 bp

All species contained:
13 protein-coding genes
22 tRNAs
2 rRNAs
Gene order was completely conserved across all species
Codon Usage Bias
Strong A/T-ending codon preference observed
Leucine codon tta showed particularly strong bias
RSCU clustering separated species into distinct codon usage groups
Natural Selection

All protein-coding genes showed:
Ka/Ks < 1
Indicates strong purifying selection

Example values:
Gene	Mean Ka/Ks
cox1	0.013
cob	0.022
atp8	0.124
Phylogenetics
High bootstrap support values (57–100)
Bayesian posterior probabilities up to 1.0
T. capricornis and T. paradussumieri consistently grouped as sister taxa

Tools & Software
Bioinformatics Tools
MitoZ
MITOS2
MEGA
RevBayes
RevGadgets
BLAST
Proksee
Programming Languages
Python
R
R Packages
ggplot2
RSCUcaller
RevGadgets
Example Outputs
Generated Analyses
Circular mitochondrial genome maps
RSCU histograms
Codon usage heatmaps
Ka/Ks bar plots
Maximum Likelihood phylogenetic trees
Bayesian phylogenetic trees
