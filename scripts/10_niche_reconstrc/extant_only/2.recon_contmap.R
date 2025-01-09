

library(ape)
library(phytools)
library(ggplot2)
library(ggpubr)



pnos_dir = "./pnos/"
  
vars = paste("pno",seq(1,12),sep="")
specieslist = read.csv("./specieslist_409.csv",header=F)$V1
names(specieslist) = read.csv("./specieslist_409.csv",header=F)$V3
specieslist

mcc = read.nexus("scripts/5_BEAST/outputs/mcctree.calib0a9f_conMO431_89k_top20_10000resample.tre")
mcc.pr = keep.tip(mcc,tip=specieslist) 
mcc.pr.la = ladderize(mcc.pr,right=T)
write.tree(mcc.pr.la, "./outputs/mcctree.calib0a9f_conMO431_89k_top20_10000resample.pr.la.tre")

mcc.pr.la = read.tree("./outputs/mcctree.calib0a9f_conMO431_89k_top20_10000resample.pr.la.tre")
# pdf("nola.pdf",height=18,width = 9)
# plot(mcc.pr)
# plot(mcc.pr.la)
# dev.off()

contMap.singleVar.MutiSamples <- function(variable=vars[1],specieslist=specieslist,value_reps=100,tree=mcc.pr.la){
  
  # (1)generate the replicates of niche values for each species 
  niche_sampled = matrix(nrow=length(specieslist),ncol=value_reps)
  rownames(niche_sampled) <- specieslist
  for (i in 1:length(specieslist)){
    pno = read.csv(paste0(pnos_dir,variable,"_",specieslist[i],".csv"),header=T)
    niche_sampled[i,] = sample(x=pno$variable, size=value_reps, replace=T, prob=pno$Model_avg)
  }
  
  # (2) Apply analysis of ancestral state estimation for a single variable
  niche_anc = matrix(nrow=tree$Nnode,ncol=value_reps)
  rownames(niche_anc) <- as.character(1:tree$Nnode+Ntip(tree))
  for (j in 1:value_reps){
    print(j)
    node_niche_size0 <- fastAnc(tree, niche_sampled[,j], vars=TRUE, CI=TRUE)
    node_niche_size <- node_niche_size0$ace
    node_niche_size
    #only retain the the ancestral value for internodes
    node_niche_size <- node_niche_size[as.character(1:tree$Nnode+Ntip(tree))]
    node_niche_size
    #restore the result of ancestral state estimation
    niche_anc[,j] = node_niche_size
  
    ## projection of the reconstruction onto the edges of the tree
    #obj <- contMap(tree, niche_sampled[,j], plot=FALSE, res=100)
    
    #Extract the internal node values for the ancestral character states
    #ii <- c(as.numeric(names(obj$tree$maps[[1]])[1]),
    #      sapply(obj$tree$maps,function(x) as.numeric(names(x)[length(x)])))
    #length(x) means the length of obj$tree$maps[[?]]
    # id = vector()
    # id2 = vector()
    # for (i in 1:length(obj$tree$maps)){
    #   id = append(id,as.numeric(names(obj$tree$maps[[i]])))
    #   id2 = append(id2,obj$tree$maps[[i]])
    # }
    # summary(id)
    # summary(id2)
    #find the trait values corresponding to each index, i.e., the node number
    #node_niche_size <- setNames(ii/(length(obj$cols)-1)*diff(obj$lims)+obj$lims[1],
    #                           c(obj$tree$edge[1,1],obj$tree$edge[,2]))
    #head(node_niche_size)
    #only retain the the ancestral value for internodes
    #node_niche_size<-node_niche_size[as.character(1:tree$Nnode+Ntip(tree))]
    #node_niche_size 
    ##Easier format to transfer to excel
    #node_niche_size<-as.matrix(node_niche_size)  
    #head(node_niche_size)
    #restore the result of ancestral state estimation
    #niche_anc[,j] = node_niche_size
  }
  result = list(niche_anc,niche_sampled)
  return(result)
}



#run ancestral state estimation for 12 representative variables
pnos_anc = list()
pnos_sampled = list()
for (i in 1:length(vars)){
  print(paste0("Conducting ancestral state estimation for pno",i))
  anc = contMap.singleVar.MutiSamples(variable=vars[i],specieslist=specieslist,value_reps=100,tree=mcc.pr.la)
  pnos_anc[[i]] = anc[[1]]
  pnos_sampled[[i]] = anc[[2]]
  print(paste0("Ancestral state estimation for pno",i," is finished!!"))
}
names(pnos_anc) <- vars
names(pnos_sampled) <- vars
save(pnos_anc, file="outputs/nicheRecon.100reps.12pnos.Rdata")
save(pnos_sampled, file="outputs/nicheSampled.100reps.12pnos.Rdata")

load("outputs/nicheRecon.100reps.12pnos.Rdata")
load("outputs/nicheSampled.100reps.12pnos.Rdata")



contMap.meanVar.plot.fan <- function(var="pno1", pnos_anc_list=pnos_anc, pnos_sampled_list=pnos_sampled, tree=mcc.pr.la){
  pno_observed_m = apply(pnos_sampled_list[[var]], 1, mean)
  anc.states_m = apply(pnos_anc_list[[var]], 1, mean)
  #fix the state for all internal nodes with the mean of estimated ancestral state values from 100 replicates
  obj  = contMap(tree, pno_observed_m, anc.states=anc.states_m,
                 plot=FALSE, res=100)
  obj<-setMap(obj,invert=TRUE) ## invert color map

  #clade labels for fan tree
  plot(obj, ftype="off",type="fan",outline=FALSE,legend=NULL,
       fsize=c(0,0.8), lwd=4)
  clade_abb = unique(names(specieslist))
  for (i in 1:length(clade_abb)){
    clade_sp = specieslist[names(specieslist)==clade_abb[i]]
    if (length(clade_sp)==1){
      print(clade_abb[i])
      arc.cladelabels(text=clade_abb[i], node=which(tree$tip.label==clade_sp),
                      ln.offset=1.02,lab.offset=1.03, orientation="horizontal", mark.node=FALSE)
    }else{
      print(clade_abb[i])
      arc.cladelabels(text=clade_abb[i], node=findMRCA(tree,clade_sp),
                      ln.offset=1.02,lab.offset=1.05, orientation="curved", mark.node=FALSE)
    }
  }
}


contMap.meanVar.plot.phylogram <- function(var="pno1", pnos_anc_list=pnos_anc, pnos_sampled_list=pnos_sampled, tree=mcc.pr.la){
  pno_observed_m = apply(pnos_sampled_list[[var]], 1, mean)
  anc.states_m = apply(pnos_anc_list[[var]], 1, mean)
  #fix the state for all internal nodes with the mean of estimated ancestral state values from 100 replicates
  obj  = contMap(tree, pno_observed_m, anc.states=anc.states_m,
                 plot=FALSE, res=100)
  obj<-setMap(obj,invert=TRUE) ## invert color map
  
  plot(obj, ftype="off",type="phylogram",outline=FALSE,legend=30,
       fsize=c(0,0.8), lwd=1.2, mar=c(3.1,0.2,0.2,0.2))
  axis(1,seq(0,90,by=10))
  title(var)
  
  #Clade labels for phylogram
  clade_abb = unique(names(specieslist))
  for (i in 1:length(clade_abb)){
    clade_sp = specieslist[names(specieslist)==clade_abb[i]]
    if (length(clade_sp)==1){
      cladelabels(text=clade_abb[i], node=which(tree$tip.label==clade_sp),
                  wing.length=0, offset=0.5, orientation="vertical")
    }else{
      cladelabels(text=clade_abb[i], node=findMRCA(tree,clade_sp),
                  wing.length=0, offset=0.5, orientation="vertical")
    }
  }
}

pdf("outputs/contMap.meanVar.100reps.pno1_12.pdf",width=8,height=15)
par(mfrow=c(2,6))
contMap.meanVar.plot.phylogram(var="pno1", pnos_anc_list=pnos_anc, pnos_sampled_list=pnos_sampled, tree=mcc.pr.la)
contMap.meanVar.plot.phylogram(var="pno2", pnos_anc_list=pnos_anc, pnos_sampled_list=pnos_sampled, tree=mcc.pr.la)
contMap.meanVar.plot.phylogram(var="pno3", pnos_anc_list=pnos_anc, pnos_sampled_list=pnos_sampled, tree=mcc.pr.la)
contMap.meanVar.plot.phylogram(var="pno4", pnos_anc_list=pnos_anc, pnos_sampled_list=pnos_sampled, tree=mcc.pr.la)
contMap.meanVar.plot.phylogram(var="pno5", pnos_anc_list=pnos_anc, pnos_sampled_list=pnos_sampled, tree=mcc.pr.la)
contMap.meanVar.plot.phylogram(var="pno6", pnos_anc_list=pnos_anc, pnos_sampled_list=pnos_sampled, tree=mcc.pr.la)
contMap.meanVar.plot.phylogram(var="pno7", pnos_anc_list=pnos_anc, pnos_sampled_list=pnos_sampled, tree=mcc.pr.la)
contMap.meanVar.plot.phylogram(var="pno8", pnos_anc_list=pnos_anc, pnos_sampled_list=pnos_sampled, tree=mcc.pr.la)
contMap.meanVar.plot.phylogram(var="pno9", pnos_anc_list=pnos_anc, pnos_sampled_list=pnos_sampled, tree=mcc.pr.la)
contMap.meanVar.plot.phylogram(var="pno10", pnos_anc_list=pnos_anc, pnos_sampled_list=pnos_sampled, tree=mcc.pr.la)
contMap.meanVar.plot.phylogram(var="pno11", pnos_anc_list=pnos_anc, pnos_sampled_list=pnos_sampled, tree=mcc.pr.la)
contMap.meanVar.plot.phylogram(var="pno12", pnos_anc_list=pnos_anc, pnos_sampled_list=pnos_sampled, tree=mcc.pr.la)
dev.off()






extract.AncValue.MRCA.PCA12 <- function(tree=mcc.pr.la,
                                        MRCA_tips=c("NLdensiflorus_K003_WA07","CEvariabilis_trans_A"),
                                        node_name="SG_Notholithocarpus",
                                        pnos_anc_list=pnos_anc){
  #extract the ancestral values of 12 variables from 100 replicates for a node of certain clade
  node_number = getMRCA(tree,tip=MRCA_tips)
  AncValue.mat = sapply(pnos_anc_list, function(x) x[as.character(node_number),])
  AncValue.df = as.data.frame(AncValue.mat)
  #perfrom PCA analysis for 12 varibles of this node for result visualization
  AncValue.pr = princomp(AncValue.df, cor=TRUE)
  print(summary(AncValue.pr, loadings=TRUE))
  AncValue.pred = predict(AncValue.pr)
  head(AncValue.pred)
  PC12 = data.frame(node=rep(node_name,nrow(AncValue.pred)), PC1=AncValue.pred[,"Comp.1"], PC2=AncValue.pred[,"Comp.2"])
  return(PC12)
}



SG_NL = extract.AncValue.MRCA.PCA12(tree=mcc.pr.la,
                                    MRCA_tips=c("NLdensiflorus_K003_WA07","CEvariabilis_trans_A"),
                                    node_name="SG_Notholithocarpus",
                                    pnos_anc_list=pnos_anc)
CG_CH = extract.AncValue.MRCA.PCA12(tree=mcc.pr.la,
                                    MRCA_tips=c("CHsempervirens_U03_WC08", "CHchrysophylla_trans"),
                                    node_name="CG_Chrysolepis",
                                    pnos_anc_list=pnos_anc)
CG_QU = extract.AncValue.MRCA.PCA12(tree=mcc.pr.la,
                                    MRCA_tips=c("LOwislizeni_F079_WF10_A", "QUpyrenaica_P009_WB11_A"),
                                    node_name="CG_subg.Quercus",
                                    pnos_anc_list=pnos_anc)
CG_CE = extract.AncValue.MRCA.PCA12(tree=mcc.pr.la,
                                    MRCA_tips=c("CEvariabilis_trans_A", "CYgaharuensis_P008_WG02_A"),
                                    node_name="CG_subg.Cerris",
                                    pnos_anc_list=pnos_anc)
CG_LP = extract.AncValue.MRCA.PCA12(tree=mcc.pr.la,
                                    MRCA_tips=c("LPcorneus_P001_WB10", "LPtaitoensis_F118_WG03"),
                                    node_name="CG_Lithocarpus",
                                    pnos_anc_list=pnos_anc)
CG_CTCP = extract.AncValue.MRCA.PCA12(tree=mcc.pr.la,
                                      MRCA_tips=c("CTsativa_F102_WA11", "CPnephelioides_F116_WH10"),
                                      node_name="CG_CTCP",
                                      pnos_anc_list=pnos_anc)

node_name="SG_Notholithocarpus"   #"NLdensiflorus_K003_WA07", "CEvariabilis_trans_A"
node_name="CG_Chrysolepis"        #"CHsempervirens_U03_WC08", "CHchrysophylla_trans"
node_name="CG_subg.Quercus"       #"LOwislizeni_F079_WF10_A", "QUpyrenaica_P009_WB11_A"
node_name="CG_subg.Cerris"        #"CEvariabilis_trans_A", "CYgaharuensis_P008_WG02_A"
node_name="CG_Lithocarpus"        #"LPcorneus_P001_WB10", "LPtaitoensis_F118_WG03"
node_name="CG_CTCP"               #"CTsativa_F102_WA11", "CPnephelioides_F116_WH10"



PCA12_NW = data.frame(node=vector(), PC1=vector(), PC2=vector())
PCA12_NW = rbind(PCA12_NW,SG_NL)
PCA12_NW = rbind(PCA12_NW,CG_CH)
PCA12_NW = rbind(PCA12_NW,CG_QU)
PCA12_NW$node = factor(PCA12_NW$node, levels=c("SG_Notholithocarpus","CG_Chrysolepis","CG_subg.Quercus"),
                          labels=c("SG_Notholithocarpus","CG_Chrysolepis","CG_subg.Quercus"))

PCA12_OW = data.frame(node=vector(), PC1=vector(), PC2=vector())
PCA12_OW = rbind(PCA12_OW,CG_CE)
PCA12_OW = rbind(PCA12_OW,CG_LP)
PCA12_OW = rbind(PCA12_OW,CG_CTCP)
PCA12_OW$node = factor(PCA12_OW$node, levels=c("CG_Lithocarpus","CG_CTCP","CG_subg.Cerris"),
                       labels=c("CG_Lithocarpus","CG_CTCP","CG_subg.Cerris"))


################################################################################################
# visualization for the PC1 and PC2 vaules for above six nodes
p1 <- ggplot(PCA12_OW, aes(x=PC1, y=PC2, color=node, alpha=node))+
  geom_point()+
  stat_ellipse(level=0.95)+
  scale_color_manual(values=c("orange","steelblue","limegreen"))+
  scale_alpha_manual(values=c(1,1,1))+
  ggtitle("a) Extant species")+
  theme(legend.position=c(0.8,0.8),
        legend.key.size=unit(10,"pt"),
        legend.title=element_blank(),
        legend.text=element_text(size=5))
  
p2 <- ggplot(PCA12_NW, aes(x=PC1, y=PC2, color=node, alpha=node))+
  geom_point()+
  stat_ellipse(level=0.95)+
  scale_color_manual(values=c("orange","steelblue","limegreen"))+
  scale_alpha_manual(values=c(1,1,1))+
  ggtitle("b) Extant species")+
  theme(legend.position=c(0.8,0.8),
        legend.key.size=unit(10,"pt"),
        legend.title=element_blank(),
        legend.text=element_text(size=5))

save(p1,p2,file="outputs/p1p2.RData")



################################################################################################


load("outputs/p1p2.RData")
pdf(file="outputs/PC12plot_12pnos_6clades.pdf",width=9,height=6.5)
ggarrange(p1,p2,nrow=2,ncol=3)
dev.off()


