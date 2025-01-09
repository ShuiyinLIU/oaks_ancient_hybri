#
#According to a mini-DB namelist which should be retained, and least percent of taxa
#in pruned gene tree, this script first filter gene trees satisfying this criterion.
#then remove tips not in mini-DB and output the pruned gene tree.
#For variety of rerooted tree, you can use "parallel" or a loop of "while" for processing 
# For example:
#while read name;do 
#Rscript ~/liushuiyin/applications/Scripts_lsy/prune_rerootedgenetrees_phylonet.R ../../namelist_miniDBplanB4_HYB.txt ~/liushuiyin/Hyb_DB/RT_ot8ot7ot2R2.V3sra.tm_431/${name}.raxml_bs.tre.rr 33 genetrees_pruned/${name}.raxml_bs.tre.rr
#done <../../genelist_cand_HYB.RT114.txt

Args = commandArgs(TRUE)
minidb_file = Args[1]
gtree_file = Args[2]
minp_taxa = Args[3]
output_file = Args[4]

#"USAGE: Rscript prune_rerootedgenetrees_phylonet.R "
#for example
#minidb_file = "Phylonetwork/namelist_miniDBplanC2_HYB.txt"
#gtree_file = "Test/gtree_rr/locus52.inclade1.ortho1.raxml_bs.tre.rr"
#minp_taxa = 33
#output_file = "Test/genetrees_pruned/locus52.inclade1.ortho1.raxml_bs.tre.rr"

library(ape)

#read namelist of mini-DB
splist_df = read.csv(minidb_file, header=F, stringsAsFactors = F)
splist = splist_df$V1
#print(splist)
#get minimum taxa number  
if ((as.numeric(minp_taxa)/100)*length(splist)>4){
  minNum_taxa <- (as.numeric(minp_taxa)/100)*length(splist)
}else{
  minNum_taxa <- 3 #i.e., at least 4 taxa should be in pruned gene tree
  }


#read rerooted gene tree
tr <- read.tree(gtree_file)
overlapped_taxa <- intersect(splist,tr$tip.label)
overlapped_taxaNum <- length(overlapped_taxa)
#filter gene tree then prune and output pruned tree
if (overlapped_taxaNum>minNum_taxa){
  tr_pruned <- keep.tip(tr,tip=overlapped_taxa)
  print(is.rooted(tr_pruned))
  write.tree(tr_pruned,file=output_file)
}else{
  print(paste(gtree_file, " has less than ",minp_taxa,"% species in miniDB!",sep=""))
}





