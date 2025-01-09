

library(raster)
library(colorRamps)



plot_md <- function(MD_dir,clade,stage,stage_fossil,if_oak, cex_pol=1, cex_mac=1.5){
  #"stage", means MD of which stage will be visualized, length of 1
  #"stage_fossil", means fossil record of which stage will be plotted on MD map, length of >=1
  print(clade)
  print(stage_fossil)
  #md <- read.asciigrid(fname=paste(MD_dir,paste(clade,"_",stage,".asc",sep=""),sep="/"))
  md <- raster(paste(MD_dir,paste(clade,"_",stage,".asc",sep=""),sep="/"))
  crs(md) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  plot(md, col=rev(blue2green2red(100)),main=paste(clade,stage, sep="_"),axes=T)
  #plot fossil record
  df <- read.csv("./FossilRecords_paleocor.csv", header=T)
  if (if_oak){
    subgen = unlist(strsplit(clade,'[.]'))[2]
    df1 <- df[(df$Genus=="Quercus")&(df$Stage_mean%in%stage_fossil),]
    df1 <- df1[df1$subgenus%in%subgen,]
    df1 <- df1[!is.na(df1$paleoLong),]
    if (nrow(df1)>=1){
      if ("macrofossil"%in%unique(df1$Type)&"pollen"%in%unique(df1$Type)){
        #macro-fossil
        fossil_m <- df1[df1$Type=="macrofossil",c('paleoLong','paleoLat')]
        coordinates(fossil_m) <- c('paleoLong','paleoLat')
        crs(fossil_m) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
        points(fossil_m, pch=3, col="black", cex=cex_mac)
        #pollen fossil
        fossil_p <- df1[df1$Type=="pollen",c('paleoLong','paleoLat')]
        coordinates(fossil_p) <- c('paleoLong','paleoLat')
        crs(fossil_p) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
        points(fossil_p, pch=16, col="black", cex=cex_pol)
      } else if ("macrofossil"%in%unique(df1$Type)&!("pollen"%in%unique(df1$Type))){
        fossil_m <- df1[df1$Type=="macrofossil",c('paleoLong','paleoLat')]
        coordinates(fossil_m) <- c('paleoLong','paleoLat')
        crs(fossil_m) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
        points(fossil_m, pch=3, col="black", cex=cex_mac)
      } else if (!("macrofossil"%in%unique(df1$Type))&("pollen"%in%unique(df1$Type))){
        fossil_p <- df1[df1$Type=="pollen",c('paleoLong','paleoLat')]
        coordinates(fossil_p) <- c('paleoLong','paleoLat')
        crs(fossil_p) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
        points(fossil_p, pch=16, col="black", cex=cex_pol)
      }
    }else{
      print("There is no fossil record for this clade and stage")
    }
  } else{
    df1 <- df[df$Genus==clade&df$Stage_mean%in%stage_fossil,]
    df1 <- df1[!is.na(df1$paleoLong),]
    if (nrow(df1)>=1){
      if ("macrofossil"%in%unique(df1$Type)&"pollen"%in%unique(df1$Type)){
        #macro-fossil
        fossil_m <- df1[df1$Type=="macrofossil",c('paleoLong','paleoLat')]
        coordinates(fossil_m) <- c('paleoLong','paleoLat')
        crs(fossil_m) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
        points(fossil_m, pch=3, col="black", cex=cex_mac)
        #pollen fossil
        fossil_p <- df1[df1$Type=="pollen",c('paleoLong','paleoLat')]
        coordinates(fossil_p) <- c('paleoLong','paleoLat')
        crs(fossil_p) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
        points(fossil_p, pch=16, col="black", cex=cex_pol)
      } else if ("macrofossil"%in%unique(df1$Type)&!("pollen"%in%unique(df1$Type))){
        fossil_m <- df1[df1$Type=="macrofossil",c('paleoLong','paleoLat')]
        coordinates(fossil_m) <- c('paleoLong','paleoLat')
        crs(fossil_m) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
        points(fossil_m, pch=3, col="black", cex=cex_mac)
      } else if (!("macrofossil"%in%unique(df1$Type))&("pollen"%in%unique(df1$Type))){
        fossil_p <- df1[df1$Type=="pollen",c('paleoLong','paleoLat')]
        coordinates(fossil_p) <- c('paleoLong','paleoLat')
        crs(fossil_p) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
        points(fossil_p, pch=16, col="black", cex=cex_pol)
      }
    }else{
      print("There is no fossil record for this clade and stage")
    }
  }
}


# set parameters
clade_list=c("Castanea","Castanopsis","Lithocarpus","subgen.Cerris","subgen.Quercus")
stage_mask_PE=c("Danian","Thanetian","Ypresian", "Lutetian", "Bartonian", "Priabonian")
stage_mask_PEO=c("Danian","Thanetian","Ypresian", "Lutetian", "Bartonian", "Priabonian", "Rupelian", "Chattian")


########## PE
#L. Paleocene, E. Eocene, M. Eocene ("Lutetian", "Bartonian"), L. Eocene
pdf("6_MD_plots/PE.L.paleocene_2L.Eocene_10w9h.pdf", width=9.8, height=9.5)
par(mfrow=c(5,4))
plot_md(MD_dir="./6_MD_raster_PE",clade=clade_list[1],stage=stage_mask_PE[2],stage_fossil=stage_mask_PE[2],if_oak=F, cex_pol=1.5, cex_mac=2)
plot_md(MD_dir="./6_MD_raster_PE",clade=clade_list[1],stage=stage_mask_PE[3],stage_fossil=stage_mask_PE[3],if_oak=F, cex_pol=1.5, cex_mac=2)
plot_md(MD_dir="./6_MD_raster_PE",clade=clade_list[1],stage=stage_mask_PE[4],stage_fossil=stage_mask_PE[4:5],if_oak=F, cex_pol=1.5, cex_mac=2)
plot_md(MD_dir="./6_MD_raster_PE",clade=clade_list[1],stage=stage_mask_PE[6],stage_fossil=stage_mask_PE[6],if_oak=F, cex_pol=1.5, cex_mac=2)

plot_md(MD_dir="./6_MD_raster_PE",clade=clade_list[2],stage=stage_mask_PE[2],stage_fossil=stage_mask_PE[2],if_oak=F, cex_pol=1.5, cex_mac=2)
plot_md(MD_dir="./6_MD_raster_PE",clade=clade_list[2],stage=stage_mask_PE[3],stage_fossil=stage_mask_PE[3],if_oak=F, cex_pol=1.5, cex_mac=2)
plot_md(MD_dir="./6_MD_raster_PE",clade=clade_list[2],stage=stage_mask_PE[4],stage_fossil=stage_mask_PE[4:5],if_oak=F, cex_pol=1.5, cex_mac=2)
plot_md(MD_dir="./6_MD_raster_PE",clade=clade_list[2],stage=stage_mask_PE[6],stage_fossil=stage_mask_PE[6],if_oak=F, cex_pol=1.5, cex_mac=2)

plot_md(MD_dir="./6_MD_raster_PE",clade=clade_list[3],stage=stage_mask_PE[2],stage_fossil=stage_mask_PE[2],if_oak=F, cex_pol=1.5, cex_mac=2)
plot_md(MD_dir="./6_MD_raster_PE",clade=clade_list[3],stage=stage_mask_PE[3],stage_fossil=stage_mask_PE[3],if_oak=F, cex_pol=1.5, cex_mac=2)
plot_md(MD_dir="./6_MD_raster_PE",clade=clade_list[3],stage=stage_mask_PE[4],stage_fossil=stage_mask_PE[4:5],if_oak=F, cex_pol=1.5, cex_mac=2)
plot_md(MD_dir="./6_MD_raster_PE",clade=clade_list[3],stage=stage_mask_PE[6],stage_fossil=stage_mask_PE[6],if_oak=F, cex_pol=1.5, cex_mac=2)

plot_md(MD_dir="./6_MD_raster_PE",clade=clade_list[4],stage=stage_mask_PE[2],stage_fossil=stage_mask_PE[2],if_oak=T, cex_pol=1.5, cex_mac=2)
plot_md(MD_dir="./6_MD_raster_PE",clade=clade_list[4],stage=stage_mask_PE[3],stage_fossil=stage_mask_PE[3],if_oak=T, cex_pol=1.5, cex_mac=2)
plot_md(MD_dir="./6_MD_raster_PE",clade=clade_list[4],stage=stage_mask_PE[4],stage_fossil=stage_mask_PE[4:5],if_oak=T, cex_pol=1.5, cex_mac=2)
plot_md(MD_dir="./6_MD_raster_PE",clade=clade_list[4],stage=stage_mask_PE[6],stage_fossil=stage_mask_PE[6],if_oak=T, cex_pol=1.5, cex_mac=2)

plot_md(MD_dir="./6_MD_raster_PE",clade=clade_list[5],stage=stage_mask_PE[2],stage_fossil=stage_mask_PE[2],if_oak=T, cex_pol=1.5, cex_mac=2)
plot_md(MD_dir="./6_MD_raster_PE",clade=clade_list[5],stage=stage_mask_PE[3],stage_fossil=stage_mask_PE[3],if_oak=T, cex_pol=1.5, cex_mac=2)
plot_md(MD_dir="./6_MD_raster_PE",clade=clade_list[5],stage=stage_mask_PE[4],stage_fossil=stage_mask_PE[4:5],if_oak=T, cex_pol=1.5, cex_mac=2)
plot_md(MD_dir="./6_MD_raster_PE",clade=clade_list[5],stage=stage_mask_PE[6],stage_fossil=stage_mask_PE[6],if_oak=T, cex_pol=1.5, cex_mac=2)
dev.off()



########## PEO
pdf("6_MD_plots/PEO.L.paleocene_2L.Oligocene_12w8h.pdf", width=12, height=8)
par(mfrow=c(5,6))
plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[1],stage=stage_mask_PEO[2],stage_fossil=stage_mask_PEO[2],if_oak=F)
plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[1],stage=stage_mask_PEO[3],stage_fossil=stage_mask_PEO[3],if_oak=F)
plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[1],stage=stage_mask_PEO[4],stage_fossil=stage_mask_PEO[4:5],if_oak=F)
plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[1],stage=stage_mask_PEO[6],stage_fossil=stage_mask_PEO[6],if_oak=F)
plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[1],stage=stage_mask_PEO[7],stage_fossil=stage_mask_PEO[7],if_oak=F)
plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[1],stage=stage_mask_PEO[8],stage_fossil=stage_mask_PEO[8],if_oak=F)

plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[2],stage=stage_mask_PEO[2],stage_fossil=stage_mask_PEO[2],if_oak=F)
plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[2],stage=stage_mask_PEO[3],stage_fossil=stage_mask_PEO[3],if_oak=F)
plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[2],stage=stage_mask_PEO[4],stage_fossil=stage_mask_PEO[4:5],if_oak=F)
plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[2],stage=stage_mask_PEO[6],stage_fossil=stage_mask_PEO[6],if_oak=F)
plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[2],stage=stage_mask_PEO[7],stage_fossil=stage_mask_PEO[7],if_oak=F)
plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[2],stage=stage_mask_PEO[8],stage_fossil=stage_mask_PEO[8],if_oak=F)

plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[3],stage=stage_mask_PEO[2],stage_fossil=stage_mask_PEO[2],if_oak=F)
plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[3],stage=stage_mask_PEO[3],stage_fossil=stage_mask_PEO[3],if_oak=F)
plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[3],stage=stage_mask_PEO[4],stage_fossil=stage_mask_PEO[4:5],if_oak=F)
plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[3],stage=stage_mask_PEO[6],stage_fossil=stage_mask_PEO[6],if_oak=F)
plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[3],stage=stage_mask_PEO[7],stage_fossil=stage_mask_PEO[7],if_oak=F)
plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[3],stage=stage_mask_PEO[8],stage_fossil=stage_mask_PEO[8],if_oak=F)

plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[4],stage=stage_mask_PEO[2],stage_fossil=stage_mask_PEO[2],if_oak=T)
plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[4],stage=stage_mask_PEO[3],stage_fossil=stage_mask_PEO[3],if_oak=T)
plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[4],stage=stage_mask_PEO[4],stage_fossil=stage_mask_PEO[4:5],if_oak=T)
plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[4],stage=stage_mask_PEO[6],stage_fossil=stage_mask_PEO[6],if_oak=T)
plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[4],stage=stage_mask_PEO[7],stage_fossil=stage_mask_PEO[7],if_oak=T)
plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[4],stage=stage_mask_PEO[8],stage_fossil=stage_mask_PEO[8],if_oak=T)

plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[5],stage=stage_mask_PEO[2],stage_fossil=stage_mask_PEO[2],if_oak=T)
plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[5],stage=stage_mask_PEO[3],stage_fossil=stage_mask_PEO[3],if_oak=T)
plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[5],stage=stage_mask_PEO[4],stage_fossil=stage_mask_PEO[4:5],if_oak=T)
plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[5],stage=stage_mask_PEO[6],stage_fossil=stage_mask_PEO[6],if_oak=T)
plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[5],stage=stage_mask_PEO[7],stage_fossil=stage_mask_PEO[7],if_oak=T)
plot_md(MD_dir="./6_MD_raster_PEO",clade=clade_list[5],stage=stage_mask_PEO[8],stage_fossil=stage_mask_PEO[8],if_oak=T)

dev.off()




