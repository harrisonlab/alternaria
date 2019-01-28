

GGtree was used to make a plot.

* Note- Tips can be found here: https://bioconnector.org/r-ggtree.html

* Note- Tips can be found here: https://bioconnector.org/r-ggtree.html

The consensus tree was downloaded to my local machine

* Note - I had to import into geneious and export again in newick format to get around polytomy branches having no branch length.
* Terminal branch lengths are meanlingless from ASTRAL and should all be set to an arbitrary value. This will be done by geneious (set to 1), but it also introduces a branch length of 2 for one isolate that needs to be corrected with sed

<!-- ```bash
cat Alt_phylogeny.consensus.scored.geneious.tre | sed 's/:2/:1/g' > Alt_phylogeny.consensus.scored.geneious2.tre
``` -->


```r
setwd("/Users/armita/Downloads/Aalt/5-gene_phylogeny")
#===============================================================================
#       Load libraries
#===============================================================================

library(ape)
library(ggplot2)

library(ggtree)
library(phangorn)
library(treeio)

tree <- read.tree("/Users/armita/Downloads/Aalt/5-gene_phylogeny/5_gene_Tree_17-07-13_consensus.newick")

mydata <- read.csv("/Users/armita/Downloads/Aalt/5-gene_phylogeny/traits.csv", stringsAsFactors=FALSE)
rownames(mydata) <- mydata$Isolate
mydata <- mydata[match(tree$tip.label,rownames(mydata)),]

t <- ggtree(tree, aes(linetype=nodes$support)) # Core tree
# Adjust terminal branch lengths:
branches <- t$data
tree$edge.length[branches$isTip] <- 1.0
#Tree <- branches$branch.length
#rescale_tree(t, branches$branch.length)

t <- t + geom_treescale(offset=-1.0, fontsize = 3) # Add scalebar
# t <- t + xlim(0, 0.025) # Add more space for labels



# Colouring labels by values in another df
t <- t %<+% mydata # Allow colouring of nodes by another df
scale_color_manual(values=c("gray39","black")) # colours as defined by col2rgb

tips <- data.frame(t$data)
tips$label <- tips$ID
t <- t + geom_tiplab(data=tips, size=3, hjust=0, offset = +0.1, align=T) +
scale_color_manual(values=c("gray39","black")) # colours as defined by col2rgb

# Add in a further set of labels
tree_mod <- data.frame(t$data)
tree_mod$label <- tips$Host
t <- t + geom_tiplab(data=tree_mod, aes(label=label), size=3, hjust=0, offset = +0.75,align=T, linetype = NULL, fontface = "italic")
# Add in a further set of labels
# tree_mod <- data.frame(t$data)
tree_mod$label2 <- tips$PCR
t <- t + geom_tiplab(data=tree_mod, aes(label=label2), size=3, hjust=0, offset = +2.25, align=T, linetype = NULL)

tips$MAT <- factor(tips$MAT)
# t <- t + geom_tippoint(data=tips, aes(shape=MAT), size=2.5, show.legend = NA)

t <- t + geom_tiplab(data=tips, aes(label=MAT), size=3, hjust=0, offset = +1.75, align=T, linetype = NULL)

# Format nodes by values
nodes <- data.frame(t$data)
#nodes <- nodes[!nodes$isTip,]
nodes$label <- as.numeric(nodes[nodes$label,])
as.numeric(nodes$label)
#nodes$label[nodes$label < 80] <- ''
nodes$support[nodes$isTip] <- 'supported'
nodes$support[(!nodes$isTip) & (nodes$label > 80)] <- 'supported'
nodes$support[(!nodes$isTip) & (nodes$label < 80)] <- 'unsupported'
nodes$support[(!nodes$isTip) & (nodes$label == '')] <- 'supported'
t <- t + aes(linetype=nodes$support)
nodes$label[nodes$label > 80] <- ''
#t <- t + geom_nodelab(data=nodes, size=2, hjust=-0.05) # colours as defined by col2rgb

# Annotate a clade with a bar line
# t <- t + geom_cladelabel(node=93, label='sect. Alternaria', align=T, colour='black', offset=1.0)
t <- t + geom_cladelabel(node=140, label='gaisen clade', align=T, colour='black', offset=3)
t <- t + geom_cladelabel(node=119, label='tenuissima clade', align=T, colour='black', offset=3)
t <- t + geom_cladelabel(node=95, label='arborescens clade', align=T, colour='black', offset=3)

t <- t + geom_cladelabel(node=92, label='', align=T, color='NA', offset=3.75, barsize = 0)

# Save as PDF and force a 'huge' size plot
t <- ggsave("redrawn.pdf", width =30, height = 40, units = "cm", limitsize = FALSE)

````
