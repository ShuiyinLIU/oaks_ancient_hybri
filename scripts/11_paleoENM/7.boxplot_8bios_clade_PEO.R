#Written at 11:23, SAT, 2022-11-26, by LIUSY in R

library(ggplot2)
library(agricolae)
library(ggpubr) #ggarrange()

# read all csv files with 8 bio-climatic values for all clades and stages
files <- list.files(path = "./2.niche_modelling_deeptime/3_8bio_extract", pattern = "csv$", full.names=T)
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

# add a coloum for "new world clade" (NW) and "old world clade" (OW)
bios_df$OWNW <-rep("OW",nrow(bios_df))
table(bios_df$Clade)
bios_df$OWNW[bios_df$Clade%in%c("Chrysolepis","Notholithocarpus","subgen.Quercus")] <- "NW"
table(bios_df$OWNW)
table(bios_df$Stage)



boxplot.bio <- function(data,bio,clade_order, medbar_size=0.3){
  # T-test for each pair of clades
  print(bio)
  df <- data[,c("Clade","OWNW",bio)]
  colnames(df)[3] <- "Niche"
  model <- aov(Niche~Clade, data=df)
  comparison <- LSD.test(model,"Clade",alpha=0.05,p.adj = c("bonferroni"),group = TRUE,main = bio)
  lsd <- data.frame(clade=rownames(comparison$groups),
                    bio_a=comparison$groups$groups,
                    #bio_m=comparison$means$Niche[match(rownames(comparison$groups),rownames(comparison$means))],
                    #bio_sd=comparison$means$std[match(rownames(comparison$groups),rownames(comparison$means))],
                    bio_Max=comparison$means$Max[match(rownames(comparison$groups),rownames(comparison$means))],
                    bio_pointNum=comparison$means$r[match(rownames(comparison$groups),rownames(comparison$means))])
  #lsd$bio_sd[is.na(lsd$bio_sd)] <- 0
  #clade as a factor
  df$Clade <- factor(df$Clade,levels = clade_order)
  lsd$clade <- factor(lsd$clade,levels = clade_order)
  # boxplot
  #windowsFonts(SH = windowsFont("Arial"))
  boxplot <- ggplot()+
    geom_boxplot(data=df, aes(x=Clade,y=Niche, fill=OWNW), size=medbar_size, outlier.shape = NA)+
    #geom_jitter(data=df, aes(x=Clade,y=Niche), width=0.3, shape=21, size=1) +
    geom_text(data=lsd, mapping=aes(x=clade,y=bio_Max+0.05*bio_Max,label=bio_a),size=3)+
    xlab("")+ylab("")+
    theme(axis.title.y = element_text(face="plain",size="4",color = "black"),
          axis.title.x = element_text(face="plain",size="4",color = "black"),
          axis.text.x =  element_text(size="6",color = "black", angle=45),
          axis.text.y =  element_text(size="6",color = "black"),
          panel.background = element_rect(fill = "transparent",color="black"),
          panel.border = element_rect(fill = "transparent",color = "black",linetype = 1),
          panel.grid =element_blank(),
          #legend.text = element_text(colour="black", size = 11, face = "plain"),
          #legend.title=element_blank(),
          legend.position = "none")
  return(boxplot)
}


#################################
## plot
stage_list <- c("Ypresian","Lutetian","Bartonian","Priabonian","Rupelian","Chattian")
output_dir <- "./7_boxplots_PEO"


#### Eocene for comparison
bios_df.E <- bios_df[bios_df$Stage%in%stage_list[1:4],]
table(bios_df.E$Clade)
colnames(bios_df.E)
clade_order=c("Castanea","Castanopsis","Lithocarpus","subgen.Cerris","subgen.Quercus","Chrysolepis")
CMM <- boxplot.bio(data=bios_df.E,bio="CMM",clade_order=clade_order)
DryMon <- boxplot.bio(data=bios_df.E,bio="DryMon",clade_order=clade_order)
MAP <- boxplot.bio(data=bios_df.E,bio="MAP",clade_order=clade_order)
MAT <- boxplot.bio(data=bios_df.E,bio="MAT",clade_order=clade_order)
WetDryMon <- boxplot.bio(data=bios_df.E,bio="WetDryMon",clade_order=clade_order)
WetMon <- boxplot.bio(data=bios_df.E,bio="WetMon",clade_order=clade_order)
WMM <- boxplot.bio(data=bios_df.E,bio="WMM",clade_order=clade_order)
WMMCMM <- boxplot.bio(data=bios_df.E,bio="WMMCMM",clade_order=clade_order)
  
pdf(paste(output_dir,"Eocene_8bios_w8h4.pdf",sep="/"),width = 8,height = 4)  
ggarrange(MAT,WMM,CMM,WMMCMM,MAP,WetMon,DryMon,WetDryMon, ncol = 4, nrow = 2)
dev.off()


#### Oligocene for comparison
bios_df.O <- bios_df[bios_df$Stage%in%stage_list[5:6],]
table(bios_df.O$Clade)
clade_order=c("Castanea","Castanopsis","Lithocarpus","subgen.Cerris","subgen.Quercus","Chrysolepis","Notholithocarpus")
CMM <- boxplot.bio(data=bios_df.O,bio="CMM",clade_order=clade_order)
DryMon <- boxplot.bio(data=bios_df.O,bio="DryMon",clade_order=clade_order)
MAP <- boxplot.bio(data=bios_df.O,bio="MAP",clade_order=clade_order)
MAT <- boxplot.bio(data=bios_df.O,bio="MAT",clade_order=clade_order)
WetDryMon <- boxplot.bio(data=bios_df.O,bio="WetDryMon",clade_order=clade_order)
WetMon <- boxplot.bio(data=bios_df.O,bio="WetMon",clade_order=clade_order)
WMM <- boxplot.bio(data=bios_df.O,bio="WMM",clade_order=clade_order)
WMMCMM <- boxplot.bio(data=bios_df.O,bio="WMMCMM",clade_order=clade_order)

pdf(paste(output_dir,"Oligocene_8bios_w8h4.pdf",sep="/"),width = 8,height = 4)  
ggarrange(MAT,WMM,CMM,WMMCMM,MAP,WetMon,DryMon,WetDryMon, ncol = 4, nrow = 2)
dev.off()


#### E.Eocene-L.Oligocene for comparison
#E. Eocene
bios_df.sub <- bios_df[bios_df$Stage%in%stage_list[1],]
table(bios_df.sub$Clade)
clade_order=c("Castanea","Castanopsis","Lithocarpus","subgen.Cerris","subgen.Quercus")
CMM1 <- boxplot.bio(data=bios_df.sub,bio="CMM",clade_order=clade_order)
DryMon1 <- boxplot.bio(data=bios_df.sub,bio="DryMon",clade_order=clade_order)
MAP1 <- boxplot.bio(data=bios_df.sub,bio="MAP",clade_order=clade_order)
MAT1 <- boxplot.bio(data=bios_df.sub,bio="MAT",clade_order=clade_order)
WetDryMon1 <- boxplot.bio(data=bios_df.sub,bio="WetDryMon",clade_order=clade_order)
WetMon1 <- boxplot.bio(data=bios_df.sub,bio="WetMon",clade_order=clade_order)
WMM1 <- boxplot.bio(data=bios_df.sub,bio="WMM",clade_order=clade_order)
WMMCMM1 <- boxplot.bio(data=bios_df.sub,bio="WMMCMM",clade_order=clade_order)

#M. Eocene
bios_df.sub <- bios_df[bios_df$Stage%in%stage_list[2:3],]
table(bios_df.sub$Clade)
clade_order=c("Castanea","Castanopsis","Lithocarpus","subgen.Cerris","subgen.Quercus","Chrysolepis")
CMM2 <- boxplot.bio(data=bios_df.sub,bio="CMM",clade_order=clade_order)
DryMon2 <- boxplot.bio(data=bios_df.sub,bio="DryMon",clade_order=clade_order)
MAP2 <- boxplot.bio(data=bios_df.sub,bio="MAP",clade_order=clade_order)
MAT2 <- boxplot.bio(data=bios_df.sub,bio="MAT",clade_order=clade_order)
WetDryMon2 <- boxplot.bio(data=bios_df.sub,bio="WetDryMon",clade_order=clade_order)
WetMon2 <- boxplot.bio(data=bios_df.sub,bio="WetMon",clade_order=clade_order)
WMM2 <- boxplot.bio(data=bios_df.sub,bio="WMM",clade_order=clade_order)
WMMCMM2 <- boxplot.bio(data=bios_df.sub,bio="WMMCMM",clade_order=clade_order)

#L. Eocene
bios_df.sub <- bios_df[bios_df$Stage%in%stage_list[4],]
table(bios_df.sub$Clade)
clade_order=c("Castanea","Castanopsis","Lithocarpus","subgen.Cerris","subgen.Quercus")
CMM3 <- boxplot.bio(data=bios_df.sub,bio="CMM",clade_order=clade_order)
DryMon3 <- boxplot.bio(data=bios_df.sub,bio="DryMon",clade_order=clade_order)
MAP3 <- boxplot.bio(data=bios_df.sub,bio="MAP",clade_order=clade_order)
MAT3 <- boxplot.bio(data=bios_df.sub,bio="MAT",clade_order=clade_order)
WetDryMon3 <- boxplot.bio(data=bios_df.sub,bio="WetDryMon",clade_order=clade_order)
WetMon3 <- boxplot.bio(data=bios_df.sub,bio="WetMon",clade_order=clade_order)
WMM3 <- boxplot.bio(data=bios_df.sub,bio="WMM",clade_order=clade_order)
WMMCMM3 <- boxplot.bio(data=bios_df.sub,bio="WMMCMM",clade_order=clade_order)

#E. Oligocene
bios_df.sub <- bios_df[bios_df$Stage%in%stage_list[5],]
table(bios_df.sub$Clade)
clade_order=c("Castanea","Castanopsis","Lithocarpus","subgen.Cerris","subgen.Quercus","Chrysolepis","Notholithocarpus")
CMM4 <- boxplot.bio(data=bios_df.sub,bio="CMM",clade_order=clade_order)
DryMon4 <- boxplot.bio(data=bios_df.sub,bio="DryMon",clade_order=clade_order)
MAP4 <- boxplot.bio(data=bios_df.sub,bio="MAP",clade_order=clade_order)
MAT4 <- boxplot.bio(data=bios_df.sub,bio="MAT",clade_order=clade_order)
WetDryMon4 <- boxplot.bio(data=bios_df.sub,bio="WetDryMon",clade_order=clade_order)
WetMon4 <- boxplot.bio(data=bios_df.sub,bio="WetMon",clade_order=clade_order)
WMM4 <- boxplot.bio(data=bios_df.sub,bio="WMM",clade_order=clade_order)
WMMCMM4 <- boxplot.bio(data=bios_df.sub,bio="WMMCMM",clade_order=clade_order)

#L. Oligocene
bios_df.sub <- bios_df[bios_df$Stage%in%stage_list[6],]
table(bios_df.sub$Clade)
clade_order=c("Castanea","Castanopsis","Lithocarpus","subgen.Cerris","subgen.Quercus")
CMM5 <- boxplot.bio(data=bios_df.sub,bio="CMM",clade_order=clade_order)
DryMon5 <- boxplot.bio(data=bios_df.sub,bio="DryMon",clade_order=clade_order)
MAP5 <- boxplot.bio(data=bios_df.sub,bio="MAP",clade_order=clade_order)
MAT5 <- boxplot.bio(data=bios_df.sub,bio="MAT",clade_order=clade_order)
WetDryMon5 <- boxplot.bio(data=bios_df.sub,bio="WetDryMon",clade_order=clade_order)
WetMon5 <- boxplot.bio(data=bios_df.sub,bio="WetMon",clade_order=clade_order)
WMM5 <- boxplot.bio(data=bios_df.sub,bio="WMM",clade_order=clade_order)
WMMCMM5 <- boxplot.bio(data=bios_df.sub,bio="WMMCMM",clade_order=clade_order)

pdf(paste(output_dir,"E.Eocene_2L.Oligocene_8bios_w8h12.pdf",sep="/"),width = 8,height = 12)  
ggarrange(MAT1,MAT2,MAT3,MAT4,MAT5,
          WMM1,WMM2,WMM3,WMM4,WMM5,
          CMM1,CMM2,CMM3,CMM4,CMM5,
          WMMCMM1,WMMCMM2,WMMCMM3,WMMCMM4,WMMCMM5,
          MAP1,MAP2,MAP3,MAP4,MAP5,
          WetMon1,WetMon2,WetMon3,WetMon4,WetMon5,
          DryMon1,DryMon2,DryMon3,DryMon4,DryMon5,
          WetDryMon1,WetDryMon2,WetDryMon3,WetDryMon4,WetDryMon5,
          ncol = 5, nrow = 8)
dev.off()




