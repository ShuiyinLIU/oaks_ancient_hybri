#####extract climate range based on fossil record#####
library(sp)
library(raster)
library(stringr)


extract.biovalue <- function(records_df, stage, NorthHem_bio_dir, output_dir, clade){
  # choose records in specific stage, base on "Stage_mean"
  records_df1 <- subset(records_df, Stage_mean==stage)
  
  if (nrow(records_df1)>=1){
    # read climate files in specific stage
    files_raw <- list.files(path = NorthHem_bio_dir, pattern = "asc$", full.names=T)
    files <- files_raw[grep(stage, files_raw)]
    print(files)
    rs <- stack(files)
    
    # check occurrences within the same grid, and randomly retain a single occurrence
    records_df1$cellNumber <- cellFromXY(rs[[1]], records_df1[,c('paleoLong','paleoLat')])
    records_df2 <- records_df1[!duplicated(records_df1$cellNumber),]
    records_df3 <- records_df2[,c('paleoLong','paleoLat')]
    sp_new <- as.data.frame(records_df3)
    coordinates(sp_new) <- c('paleoLong','paleoLat')
    
    # extract bio climate and write it to your direction
    #If df=TRUE, the results are returned as a 'dataframe'##?CRS is NA
    bio_extract <- extract(rs, sp_new, method='simple', buffer=1000, fun=mean, df=TRUE)
    bio_extract <- na.omit(bio_extract)
    write.csv(bio_extract, file = paste(output_dir,paste(paste(clade,stage,sep="_"), "8bio_extract.csv",sep="_"), sep="/"),row.names=FALSE)
    
  }else{
    print(paste("This clade has no fossil record in stage", stage))
  }
}


#read fossil records 
sp <- read.csv("FossilRecords_paleocor.csv", header=TRUE, sep=",")
colnames(sp)
unique(sp$Genus)



stage_list <- c("Thanetian","Ypresian","Lutetian", "Bartonian","Priabonian", "Rupelian","Chattian")
genus_list <- c("Castanea", "Castanopsis", "Lithocarpus", "Chrysolepis", "Notholithocarpus")

## genus
for (i in 1:length(genus_list)){
  print(genus_list[i])
  sp1 <- subset(sp, Genus==genus_list[i])
  for (j in 1:length(stage_list)){
    print(stage_list[j])
    extract.biovalue(records_df=sp1, stage=stage_list[j], NorthHem_bio_dir="./1_rasterdata_NorthHem_asc_5min",
                     output_dir="./3_8bio_extract", clade=genus_list[i])
  }
}


## subgenus of oaks
sp_oak <- subset(sp, Genus=="Quercus")
sp_oak_cerris <- subset(sp_oak, subgenus=="Cerris")
sp_oak_quercus <- subset(sp_oak, subgenus=="Quercus")
for (j in 1:length(stage_list)){
  print(stage_list[j])
  extract.biovalue(records_df=sp_oak_cerris, stage=stage_list[j], NorthHem_bio_dir="./1_rasterdata_NorthHem_asc_5min",
                   output_dir="./3_8bio_extract", clade="subgen.Cerris")
}

for (j in 1:length(stage_list)){
  print(stage_list[j])
  extract.biovalue(records_df=sp_oak_quercus, stage=stage_list[j], NorthHem_bio_dir="./1_rasterdata_NorthHem_asc_5min",
                   output_dir="./3_8bio_extract", clade="subgen.Quercus")
}



########################################################
# two coordinates of Castanea and Castanopsis in late Paleocene ("Thanetian") all with NAs from teydm1 (60.5Ma) paleoclimates
# below, we change to extract values from teydn1 (66.0Ma) paleoclimates

files <- list.files(path = "H:/teyd_mask_asc_5min/teydn1", pattern = "asc$", full.names=T)
files <- files[c(1,3,5,6,8,10,12,13)]
files
rs <- stack(files)

# Castanea
records_df <- subset(sp, Genus=="Castanea")
records_df1 <- subset(records_df, Stage_mean=="Thanetian")
records_df1$cellNumber <- cellFromXY(rs[[1]], records_df1[,c('paleoLong','paleoLat')])
records_df2 <- records_df1[!duplicated(records_df1$cellNumber),]
records_df3 <- records_df2[,c('paleoLong','paleoLat')]
sp_new <- as.data.frame(records_df3)
coordinates(sp_new) <- c('paleoLong','paleoLat')
bio_extract <- extract(rs, sp_new, method='simple', buffer=1000, fun=mean, df=TRUE)
bio_extract <- na.omit(bio_extract)
write.csv(bio_extract, file = "./3_8bio_extract/Castanea_Danian_8bio_extract.csv",row.names=FALSE)
# one with NAs, another with values

# Castanopsis
records_df <- subset(sp, Genus=="Castanopsis")
records_df1 <- subset(records_df, Stage_mean=="Thanetian")
records_df1$cellNumber <- cellFromXY(rs[[1]], records_df1[,c('paleoLong','paleoLat')])
records_df2 <- records_df1[!duplicated(records_df1$cellNumber),]
records_df3 <- records_df2[,c('paleoLong','paleoLat')]
sp_new <- as.data.frame(records_df3)
coordinates(sp_new) <- c('paleoLong','paleoLat')
bio_extract <- extract(rs, sp_new, method='simple', buffer=1000, fun=mean, df=TRUE)
bio_extract <- na.omit(bio_extract)
write.csv(bio_extract, file = "./3_8bio_extract/Castanopsis_Danian_8bio_extract.csv",row.names=FALSE)
# two with NAs.

