# Attempting Rev Bayes for BNFO620

#Installing
install.packages("devtools")
devtools::install_github("cmt2/RevGadgets")


#packages
library(RevGadgets)
library(coda)
library(ggplot2)
library(ggtree)
library(grid)
library(gridExtra)
library(treeio)
library(ape)


tree <- readTrees("/Users/sundaywright/output/concatenated_MAP.tre")

#outgroup 
tree_rooted <- rerootPhylo(tree = tree, outgroup = "G_borealis")

#
tree_rooted[[1]][[1]]@phylo$tip.label <- gsub("T_Capricornis", "T_capricornis", tree_rooted[[1]][[1]]@phylo$tip.label)

# create the plot of the rooted tree
plot <- plotTree(tree = tree_rooted,
                 # make tree lines more narrow
                 line_width = 0.5,
                 # italicize tip labels 
                 tip_labels_italics = TRUE,
                 
                 node_labels = 'posterior',
                 node_labels_size = 2.75,
                 node_labels_offset = 0.005,
                 node_labels_digits = 3
                 )

# add scale bar to the tree and plot with ggtree
plot <- plot + theme_tree() + geom_treescale(x = 0, y = 0.5, width = 0.1, offset = 0.5) +   theme(
  axis.line.x = element_blank(),    #removes line
  axis.text.x = element_blank(),    #removes text
  axis.ticks.x = element_blank(),   #removes ticks
  axis.title.x = element_blank()    #removes title
)

suppressWarnings(print(plot))