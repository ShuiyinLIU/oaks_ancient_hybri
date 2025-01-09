
library(ggplot2)
library(ggpubr)

################################################################################################
#visualization for the PC1 and PC2 of niche of Eocene fossils for several focused clades


# read all csv files with 8 bio-climatic values for all clades and stages
files <- list.files(path = "11_paleoENM/4_8bio_extract", pattern = "csv$", full.names=T)
bios_df <- as.data.frame(matrix(nrow=0,ncol=11))
names <- c("ID","CMM","DryMon","MAP","MAT","WetDryMon","WetMon","WMM","WMMCMM","Clade","Stage")
colnames(bios_df) <- names 

for (i in 1:length(files)){
  print(files[i])
  df <- read.csv(files[i], header=T)
  clade = head(unlist(strsplit(tail(unlist(strsplit(files[i],"/")),1),"_")),1)
  stage = tail(head(unlist(strsplit(tail(unlist(strsplit(files[i],"/")),1),"_")),2),1)
  df$clade <- rep(clade,nrow(df))
  df$stage <- rep(stage,nrow(df))
  colnames(df) <- names
  bios_df <- rbind(bios_df,df)
}


stage_list <- c("Ypresian","Lutetian","Bartonian","Priabonian","Rupelian","Chattian")


#### Eocene for comparison
bios_df.E <- bios_df[bios_df$Stage%in%stage_list[1:4],]
table(bios_df.E$Clade)


fossilNiche.PCA12 <- function(niche_df=bios_df.E, clade="Castanea"){
  #extract the niche values of 8 variables for all fossil record of certain clade
  niche_df.pr = niche_df[niche_df$Clade==clade,2:9]
  #perform PCA analysis for 8 variables of this clade for result visualization
  niche.pr = princomp(niche_df.pr, cor=TRUE)
  print(summary(niche.pr, loadings=TRUE))
  niche.pred = predict(niche.pr)
  head(niche.pred)
  PC12 = data.frame(node=rep(clade,nrow(niche.pred)), PC1=niche.pred[,"Comp.1"], PC2=niche.pred[,"Comp.2"])
  return(PC12)
}


CT = fossilNiche.PCA12(niche_df=bios_df.E, clade="Castanea")
CP = fossilNiche.PCA12(niche_df=bios_df.E, clade="Castanopsis")
LP = fossilNiche.PCA12(niche_df=bios_df.E, clade="Lithocarpus")
CE = fossilNiche.PCA12(niche_df=bios_df.E, clade="subgen.Cerris")
QU = fossilNiche.PCA12(niche_df=bios_df.E, clade="subgen.Quercus")


PCA12_NW = data.frame(node=vector(), PC1=vector(), PC2=vector())
PCA12_NW = rbind(PCA12_NW,QU)


PCA12_OW = data.frame(node=vector(), PC1=vector(), PC2=vector())
PCA12_OW = rbind(PCA12_OW,CE)
PCA12_OW = rbind(PCA12_OW,LP)
PCA12_OW = rbind(PCA12_OW,CP)
PCA12_OW = rbind(PCA12_OW,CT)
PCA12_OW$node = factor(PCA12_OW$node, levels=c("Lithocarpus","Castanea","Castanopsis","subgen.Cerris"),
                       labels=c("Lithocarpus","Castanea","Castanopsis","subgen.Cerris"))


# visualization for the PC1 and PC2 vaules for above 5 clades
p5 <- ggplot(PCA12_OW, aes(x=PC1, y=PC2, color=node, alpha=node))+
  geom_point()+
  stat_ellipse(level=0.95)+
  scale_color_manual(values=c("orange","steelblue","gray","limegreen"))+
  scale_alpha_manual(values=c(1,1,1,1))+
  ggtitle("e) Eocene fossils")+
  theme(legend.position=c(0.8,0.8),
        legend.key.size=unit(10,"pt"),
        legend.title=element_blank(),
        legend.text=element_text(size=5))

p6 <- ggplot(PCA12_NW, aes(x=PC1, y=PC2, color=node, alpha=node))+
  geom_point()+
  stat_ellipse(level=0.95)+
  scale_color_manual(values=c("limegreen"))+
  scale_alpha_manual(values=c(1))+
  ggtitle("f) Eocene fossils")+
  theme(legend.position=c(0.8,0.8),
        legend.key.size=unit(10,"pt"),
        legend.title=element_blank(),
        legend.text=element_text(size=5))

save(p5,p6,file="./p5p6.RData")


load("./extant_only/outputs/p1p2.RData")
load("./p5p6.RData")
pdf(file="PC12plot_12pnos_6clades.pdf",width=9,height=6.5)
ggarrange(p1,p2,p5,p6,nrow=2,ncol=3)
dev.off()
