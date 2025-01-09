#created at 09:22 2023-08-09 by LIUSY in R

library(ape)
library(phytools)
library(ggplot2)
library(ggpubr)



list.files("outputs/")
df_ex = read.csv("outputs/geo_df4_climatesCurrent.csv", header=T)
df_fo = read.csv("outputs/Fossil_paleocor_climates.csv", header=T)
#check the max and min values for each variable
summary(df_ex)
summary(df_fo)



################################################################################
################################################################################
### generate a data frame to contain climate values of 8 variable of 409 extant
### and 66 fossil species 

df_specieslist <- read.csv("./FBD_specieslist_409_66.csv",header=T)
specieslist <- df_specieslist$renameSequencingData[df_specieslist$nicheRecon=="Y"]
names(specieslist) = df_specieslist$gen_sect.abb[df_specieslist$nicheRecon=="Y"]
specieslist

### (1) read and handle dated FBD tree
# mcc = read.nexus("mcc_FBD_agerange_9c6m6_1b1.7combined.tre")
# mcc.pr = keep.tip(mcc,tip=df_specieslist$FBD_tiplabel[df_specieslist$nicheRecon=="Y"]) 
# mcc.pr.la = ladderize(mcc.pr,right=T)
# 
# tip_re <- df_specieslist$renameSequencingData[df_specieslist$nicheRecon=="Y"]
# names(tip_re) = df_specieslist$FBD_tiplabel[df_specieslist$nicheRecon=="Y"]
# mcc.pr.la$tip.label <- tip_re[mcc.pr.la$tip.label]
# write.tree(mcc.pr.la, "./outputs/mcc_FBD_agerange_9c6m6_1b1.7combined.pr.la.newick")

tre = read.tree("./outputs/mcc_FBD_agerange_9c6m6_1b1.7combined.pr.la.newick")


### (2) generate a data frame to combine values for extant and fossil species
vars <- c("CMM","DryMon","MAP","MAT","WetDryMon","WetMon","WMM","WMMCMM")

# check occurrences within the same grid, and randomly retain a single occurrence
colnames(df_ex)
df_ex_sub <- df_ex[df_ex$dup.species_cell==F,c("renameSequencingData",vars)]

colnames(df_fo)
df_fo$Bioregion_singleFossil[is.na(df_fo$Bioregion_singleFossil)] <- "N"
df_fo_sub <- df_fo[df_fo$dup.fossil_time_cell==F,c("NicheRecon_tip",vars)]
df_fo_sub_oldest <- df_fo[df_fo$dup.fossil_time_cell==F&df_fo$Bioregion_singleFossil!="",c("NicheRecon_tip",vars)]
colnames(df_fo_sub) <- colnames(df_fo_sub_oldest) <- c("renameSequencingData",vars)

data <- rbind(df_ex_sub,df_fo_sub)
data.oldest <- rbind(df_ex_sub,df_fo_sub_oldest)



################################################################################
################################################################################
### write functions for 100 randomly sampling of values and niche reconstruction

contMap.singleVar.MutiSamples <- function(data_value=data,variable=vars[1],specieslist=specieslist,value_reps=100,tree=tre){
  
  # (1)generate the replicates of niche values for each species 
  niche_sampled = matrix(nrow=length(specieslist),ncol=value_reps)
  rownames(niche_sampled) <- specieslist
  for (i in 1:length(specieslist)){
    pno = data_value[data_value$renameSequencingData==specieslist[i],variable]
    niche_sampled[i,] = sample(x=pno, size=value_reps, replace=T)
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
  }
  result = list(niche_anc,niche_sampled)
  return(result)
}



contMap.meanVar.plot.phylogram <- function(var=vars[1], pnos_anc_list=pnos_anc, pnos_sampled_list=pnos_sampled, tree=tre){
  pno_observed_m = apply(pnos_sampled_list[[var]], 1, mean)
  anc.states_m = apply(pnos_anc_list[[var]], 1, mean)
  #fix the state for all internal nodes with the mean of estimated ancestral state values from 100 replicates
  obj  = contMap(tree, pno_observed_m, anc.states=anc.states_m,
                 plot=FALSE, res=100)
  obj <- setMap(obj,invert=TRUE) ## invert color map
  
  plot(obj, ftype="off",type="phylogram",outline=FALSE,legend=50,
       fsize=c(0,0.8), lwd=1, mar=c(3.1,0.2,0.2,1.5))
  axis(1,seq(0,100,by=10))
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



extract.AncValue.MRCA.PCA12 <- function(tree=tre,
                                        MRCA_tips=c("NLdensiflorus_K003_WA07","CEvariabilis_trans_A"),
                                        node_name="SG_Notholithocarpus",
                                        pnos_anc_list=pnos_anc){
  #extract the ancestral values of 8 variables from 100 replicates for a node of certain clade
  node_number = getMRCA(tree,tip=MRCA_tips)
  AncValue.mat = sapply(pnos_anc_list, function(x) x[as.character(node_number),])
  AncValue.df = as.data.frame(AncValue.mat)
  #perfrom PCA analysis for 8 variables of this node for result visualization
  AncValue.pr = princomp(AncValue.df, cor=TRUE)
  print(summary(AncValue.pr, loadings=TRUE))
  AncValue.pred = predict(AncValue.pr)
  head(AncValue.pred)
  PC12 = data.frame(node=rep(node_name,nrow(AncValue.pred)), PC1=AncValue.pred[,"Comp.1"], PC2=AncValue.pred[,"Comp.2"])
  return(PC12)
}




################################################################################
################################################################################
### run ancestral state estimation for 8 representative variables

### use the climate value of all records for 66 fossil species
pnos_anc = list()
pnos_sampled = list()
for (i in 1:length(vars)){
  print(paste0("Conducting ancestral state estimation for ",i))
  anc = contMap.singleVar.MutiSamples(data_value=data,variable=vars[i],specieslist=specieslist,value_reps=100,tree=tre)
  pnos_anc[[i]] = anc[[1]]
  pnos_sampled[[i]] = anc[[2]]
  print(paste0("Ancestral state estimation for ",i," is finished!!"))
}
names(pnos_anc) <- names(pnos_sampled) <- vars
save(pnos_anc, file="outputs/nicheRecon.100reps.8pnos.Rdata")
save(pnos_sampled, file="outputs/nicheSampled.100reps.8pnos.Rdata")

load("outputs/nicheRecon.100reps.8pnos.Rdata")
load("outputs/nicheSampled.100reps.8pnos.Rdata")


### use the climate value of the oldest record for 66 fossil species
pnos_anc.oldest = list()
pnos_sampled.oldest = list()
for (i in 1:length(vars)){
  print(paste0("Conducting ancestral state estimation for ",i))
  anc = contMap.singleVar.MutiSamples(data_value=data.oldest,variable=vars[i],specieslist=specieslist,value_reps=100,tree=tre)
  pnos_anc.oldest[[i]] = anc[[1]]
  pnos_sampled.oldest[[i]] = anc[[2]]
  print(paste0("Ancestral state estimation for ",i," is finished!!"))
}
names(pnos_anc.oldest) <- names(pnos_sampled.oldest) <- vars
save(pnos_anc.oldest, file="outputs/nicheRecon.100reps.8pnos.oldest.Rdata")
save(pnos_sampled.oldest, file="outputs/nicheSampled.100reps.8pnos.oldest.Rdata")

load("outputs/nicheRecon.100reps.8pnos.oldest.Rdata")
load("outputs/nicheSampled.100reps.8pnos.oldest.Rdata")



################################################################################
################################################################################
### plot the estimated mean ancestral value for each variable

### use the climate value of all records for 66 fossil species
pdffn = "outputs/contMap.meanVar.100reps.pno1_8.pdf"
pdf(pdffn,width=8,height=15)
par(mfrow=c(2,4))
for (i in 1:length(vars)){
  print(i)
  contMap.meanVar.plot.phylogram(var=vars[i], pnos_anc_list=pnos_anc, pnos_sampled_list=pnos_sampled, tree=tre)
}

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it


### use the climate value of the oldest record for 66 fossil species
pdffn = "outputs/contMap.meanVar.100reps.pno1_8.oldest.pdf"
pdf(pdffn,width=8,height=15)
par(mfrow=c(2,4))
for (i in 1:length(vars)){
  contMap.meanVar.plot.phylogram(var=vars[i], pnos_anc_list=pnos_anc.oldest, pnos_sampled_list=pnos_sampled.oldest, tree=tre)
}

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it



################################################################################
################################################################################
### plot ancestral niche for main clades at PC1-PC2 niche-space

plot.PC12.Mainnode <- function(pnos_anc=pnos_anc){
  SG_NL = extract.AncValue.MRCA.PCA12(tree=tre,
                                      MRCA_tips=c("NLdensiflorus_K003_WA07","CEvariabilis_trans_A"),
                                      node_name="SG_Notholithocarpus",
                                      pnos_anc_list=pnos_anc)
  CG_CH = extract.AncValue.MRCA.PCA12(tree=tre,
                                      MRCA_tips=c("CHsempervirens_U03_WC08", "CHchrysophylla_trans"),
                                      node_name="CG_Chrysolepis",
                                      pnos_anc_list=pnos_anc)
  CG_QU = extract.AncValue.MRCA.PCA12(tree=tre,
                                      MRCA_tips=c("LOwislizeni_F079_WF10_A", "QUpyrenaica_P009_WB11_A"),
                                      node_name="CG_subg.Quercus",
                                      pnos_anc_list=pnos_anc)
  CG_CE = extract.AncValue.MRCA.PCA12(tree=tre,
                                      MRCA_tips=c("CEvariabilis_trans_A", "CYgaharuensis_P008_WG02_A"),
                                      node_name="CG_subg.Cerris",
                                      pnos_anc_list=pnos_anc)
  CG_LP = extract.AncValue.MRCA.PCA12(tree=tre,
                                      MRCA_tips=c("LPcorneus_P001_WB10", "LPtaitoensis_F118_WG03"),
                                      node_name="CG_Lithocarpus",
                                      pnos_anc_list=pnos_anc)
  CG_CTCP = extract.AncValue.MRCA.PCA12(tree=tre,
                                        MRCA_tips=c("CTsativa_F102_WA11", "CPnephelioides_F116_WH10"),
                                        node_name="CG_CTCP",
                                        pnos_anc_list=pnos_anc)
  PCA12_NW = rbind(SG_NL,CG_CH)
  PCA12_NW = rbind(PCA12_NW,CG_QU)
  PCA12_NW$node = factor(PCA12_NW$node, levels=c("SG_Notholithocarpus","CG_Chrysolepis","CG_subg.Quercus"),
                         labels=c("SG_Notholithocarpus","CG_Chrysolepis","CG_subg.Quercus"))
  
  PCA12_OW = rbind(CG_CE,CG_LP)
  PCA12_OW = rbind(PCA12_OW,CG_CTCP)
  PCA12_OW$node = factor(PCA12_OW$node, levels=c("CG_Lithocarpus","CG_CTCP","CG_subg.Cerris"),
                         labels=c("CG_Lithocarpus","CG_CTCP","CG_subg.Cerris"))
  result = list(PCA12_NW,PCA12_OW)
  return(result)
}


################################################################################################
# visualization for the PC1 and PC2 vaules for above six nodes

PC12 <- plot.PC12.Mainnode(pnos_anc=pnos_anc)
  
p3 <- ggplot(PC12[[2]], aes(x=PC1, y=PC2, color=node, alpha=node))+
  geom_point()+
  stat_ellipse(level=0.95)+
  scale_color_manual(values=c("orange","steelblue","limegreen"))+
  scale_alpha_manual(values=c(1,1,1))+
  ggtitle("c) Extant+fossil species")+
  theme(legend.position=c(0.8,0.8),
        legend.key.size=unit(10,"pt"),
        legend.title=element_blank(),
        legend.text=element_text(size=5))
  
p4 <- ggplot(PC12[[1]], aes(x=PC1, y=PC2, color=node, alpha=node))+
  geom_point()+
  stat_ellipse(level=0.95)+
  scale_color_manual(values=c("orange","steelblue","limegreen"))+
  scale_alpha_manual(values=c(1,1,1))+
  ggtitle("d) Extant+fossil species")+
  theme(legend.position=c(0.8,0.8),
        legend.key.size=unit(10,"pt"),
        legend.title=element_blank(),
        legend.text=element_text(size=5))

save(p3,p4,file="outputs/p3p4.RData")


PC12 <- plot.PC12.Mainnode(pnos_anc=pnos_anc.oldest)

p3 <- ggplot(PC12[[2]], aes(x=PC1, y=PC2, color=node, alpha=node))+
  geom_point()+
  stat_ellipse(level=0.95)+
  scale_color_manual(values=c("orange","steelblue","limegreen"))+
  scale_alpha_manual(values=c(1,1,1))+
  ggtitle("c) Extant+fossil species")+
  theme(legend.position=c(0.8,0.8),
        legend.key.size=unit(10,"pt"),
        legend.title=element_blank(),
        legend.text=element_text(size=5))

p4 <- ggplot(PC12[[1]], aes(x=PC1, y=PC2, color=node, alpha=node))+
  geom_point()+
  stat_ellipse(level=0.95)+
  scale_color_manual(values=c("orange","steelblue","limegreen"))+
  scale_alpha_manual(values=c(1,1,1))+
  ggtitle("d) Extant+fossil species")+
  theme(legend.position=c(0.8,0.8),
        legend.key.size=unit(10,"pt"),
        legend.title=element_blank(),
        legend.text=element_text(size=5))

save(p3,p4,file="outputs/p3p4.oldest.RData")



################################################################################################
#visualization for the PC1 and PC2 of niche of Eocene fossils for several focused clades

load("outputs/p3p4.RData")
pdf(file="outputs/PC12plot_8pnos_6clades.pdf",width=9,height=6.5)
ggarrange(p3,p4,nrow=2,ncol=3)
dev.off()

load("outputs/p3p4.oldest.RData")
pdf(file="outputs/PC12plot_8pnos_6clades.oldest.pdf",width=9,height=6.5)
ggarrange(p3,p4,nrow=2,ncol=3)
dev.off()

