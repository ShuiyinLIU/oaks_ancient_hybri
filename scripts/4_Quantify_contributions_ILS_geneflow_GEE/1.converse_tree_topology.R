
#USAGE: Rscript converse_tree_topology.R <treefile> <outputfile>

library(ape)

Args <- commandArgs(TRUE)
treefile <- Args[1]
outputfile <- Args[2]

#treefile <- "RAxML_bipartitions.conRT89_MOot.V1_2821.tre.rr"
#outputfile <- "topology_RAxML_bipartitions.conRT89_MOot.V1_2821.tre.rr"
tr <- read.tree(treefile)

tr_topology <- tr

tr_topology$edge.length <- tr_topology$node.label <- NULL

write.tree(tr_topology,outputfile)
