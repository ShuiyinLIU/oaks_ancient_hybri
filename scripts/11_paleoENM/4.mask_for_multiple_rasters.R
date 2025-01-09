#####extract the common geographical range of 8 bio-climate variables#####
library(raster)



get.shared.areas <- function(biovalue_dir, clade, stage_included, NorthHem_bio_dir, stage_mask, output_dir){
  #"stage_included", means the climate range of which epochs will be used to extract max-min value for ENM
  #"stage_mask", means which epochs will be masked with the max-min value for ENM
  # read climate range
  colnames <- list("ID","cmm","drymon","map","mat","wetdrymon","wetmon","wmm","wmmcmm")
  greenhouse <- as.data.frame(matrix(nrow=0,ncol=9))
  colnames(greenhouse)<-colnames
  for (i in 1:length(stage_included)){
    print(stage_included[i])
    if (stage_included[i]=="Danian"){
      climate_df <- read.csv(paste(biovalue_dir,paste0(clade,"_Danian_8bio_extract.csv"),sep="/"));colnames(climate_df)<-colnames
    } else if (stage_included[i]=="Thanetian"){
      climate_df <- read.csv(paste(biovalue_dir,paste0(clade,"_Thanetian_8bio_extract.csv"),sep="/"));colnames(climate_df)<-colnames
    } else if (stage_included[i]=="Ypresian"){
      climate_df <- read.csv(paste(biovalue_dir,paste0(clade,"_Ypresian_8bio_extract.csv"),sep="/"));colnames(climate_df)<-colnames
    } else if (stage_included[i]=="Lutetian"){
      climate_df <- read.csv(paste(biovalue_dir,paste0(clade,"_Lutetian_8bio_extract.csv"),sep="/"));colnames(climate_df)<-colnames
    } else if (stage_included[i]=="Bartonian"){
      climate_df <- read.csv(paste(biovalue_dir,paste0(clade,"_Bartonian_8bio_extract.csv"),sep="/"));colnames(climate_df)<-colnames
    } else if (stage_included[i]=="Priabonian"){
      climate_df <- read.csv(paste(biovalue_dir,paste0(clade,"_Priabonian_8bio_extract.csv"),sep="/"));colnames(climate_df)<-colnames
    } else if (stage_included[i]=="Rupelian"){
      climate_df <- read.csv(paste(biovalue_dir,paste0(clade,"_Rupelian_8bio_extract.csv"),sep="/"));colnames(climate_df)<-colnames
    } else if (stage_included[i]=="Chattian"){
      climate_df <- read.csv(paste(biovalue_dir,paste0(clade,"_Chattian_8bio_extract.csv"),sep="/"));colnames(climate_df)<-colnames
    }
    greenhouse <- rbind(greenhouse, climate_df)
  }

  # extract the max and min value of each bio-climate variable
  max <- max(greenhouse[,2])
  min <- min(greenhouse[,2])
  for (j in 2:9) {
    max[j-1] <- max(greenhouse[,j])
    min[j-1] <- min(greenhouse[,j])
  }
  print(max)
  print(min)
  
  # read climate files to a stack
  # change value out of climate range to NA, so only keep values within the range
  files_raw <- list.files(path = NorthHem_bio_dir, pattern = "asc$", full.names=T)
  for (h in 1:length(stage_mask)){
    # (1)read for a specific stage
    files <- files_raw[grep(stage_mask[h], files_raw)]
    print(files)
    rs <- stack(files)
    # (2)mask
    for (m in 1:length(rs@layers)) {
      rs[[m]][rs[[m]]<min[m]]<-NA
      rs[[m]][rs[[m]]>max[m]]<-NA
    }
    print("Mask bio layers with NA finished!")
    # (3)mask/intersect, choose the common range in all climate data
    bio_m <- mask(rs[[1]],rs[[2]])
    #plot(bio_m)
    for (n in 3:length(rs@layers)) {
      bio_m <- mask(bio_m,rs[[n]])
    }
    print("Intersect bio layers with each other finished!")
    # (4)extract 8 bio climate variables base on the mask range
    bio <- paste(clade, names(rs), sep="_")
    for (o in 1:8){
      bio_end <- mask(rs[[o]],bio_m)
      crs(bio_end) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
      writeRaster(bio_end, paste(output_dir,bio[o],sep="/"), format="ascii", overwrite=T)
    }
    print("Save 8 bio layers with mask range finished!")
  }
}



stage_mask_PE=c("Danian","Thanetian","Ypresian", "Lutetian", "Bartonian", "Priabonian")
stage_mask_PEO=c("Danian","Thanetian","Ypresian", "Lutetian", "Bartonian", "Priabonian", "Rupelian", "Chattian")

## (A) use max-min climate value of Paleocene+Eocene(PE) to mask all stages of Paleocene and Eocene layers
# only "Castanea" has "Danian" climate value.
# "Chrysolepis" (2 records) and "Notholithocarpus" (1 record) have less than three occurrences with values, so will not be
# masked and calculated MD-distance for ENM.
get.shared.areas(biovalue_dir="./3_8bio_extract", clade="Castanea",
                 stage_included=c("Danian","Ypresian", "Lutetian", "Bartonian", "Priabonian"),
                 NorthHem_bio_dir="./1_rasterdata_NorthHem_asc_5min",
                 stage_mask=stage_mask_PE, output_dir="./4_masked_bio_PE")
get.shared.areas(biovalue_dir="./3_8bio_extract", clade="Castanopsis",
                 stage_included=c("Ypresian", "Lutetian", "Bartonian", "Priabonian"),
                 NorthHem_bio_dir="./1_rasterdata_NorthHem_asc_5min",
                 stage_mask=stage_mask_PE, output_dir="./4_masked_bio_PE")
# "Lithocarpus" have no record in "Bartonian"
get.shared.areas(biovalue_dir="./3_8bio_extract", clade="Lithocarpus",
                 stage_included=c("Ypresian", "Lutetian", "Priabonian"),
                 NorthHem_bio_dir="./1_rasterdata_NorthHem_asc_5min",
                 stage_mask=stage_mask_PE, output_dir="./4_masked_bio_PE")
get.shared.areas(biovalue_dir="./3_8bio_extract", clade="subgen.Cerris",
                 stage_included=c("Ypresian", "Lutetian", "Bartonian", "Priabonian"),
                 NorthHem_bio_dir="./1_rasterdata_NorthHem_asc_5min",
                 stage_mask=stage_mask_PE, output_dir="./4_masked_bio_PE")
get.shared.areas(biovalue_dir="./3_8bio_extract", clade="subgen.Quercus",
                 stage_included=c("Ypresian", "Lutetian", "Bartonian", "Priabonian"),
                 NorthHem_bio_dir="./1_rasterdata_NorthHem_asc_5min",
                 stage_mask=stage_mask_PE, output_dir="./4_masked_bio_PE")

## (B) use max-min climate value of Paleocene+Eocene_Oligocene(PEO) to mask all stages of Paleocene, Eocene, and Oligocene layers
get.shared.areas(biovalue_dir="./3_8bio_extract", clade="Castanea",
                 stage_included=c("Danian","Ypresian", "Lutetian", "Bartonian", "Priabonian", "Rupelian", "Chattian"),
                 NorthHem_bio_dir="./1_rasterdata_NorthHem_asc_5min",
                 stage_mask=stage_mask_PEO, output_dir="./4_masked_bio_PEO")
get.shared.areas(biovalue_dir="./3_8bio_extract", clade="Castanopsis",
                 stage_included=c("Ypresian", "Lutetian", "Bartonian", "Priabonian", "Rupelian", "Chattian"),
                 NorthHem_bio_dir="./1_rasterdata_NorthHem_asc_5min",
                 stage_mask=stage_mask_PEO, output_dir="./4_masked_bio_PEO")
get.shared.areas(biovalue_dir="./3_8bio_extract", clade="Lithocarpus",
                 stage_included=c("Ypresian", "Lutetian", "Priabonian", "Rupelian", "Chattian"),
                 NorthHem_bio_dir="./1_rasterdata_NorthHem_asc_5min",
                 stage_mask=stage_mask_PEO, output_dir="./4_masked_bio_PEO")
get.shared.areas(biovalue_dir="./3_8bio_extract", clade="subgen.Cerris",
                 stage_included=c("Ypresian", "Lutetian", "Bartonian", "Priabonian", "Rupelian", "Chattian"),
                 NorthHem_bio_dir="./1_rasterdata_NorthHem_asc_5min",
                 stage_mask=stage_mask_PEO, output_dir="./4_masked_bio_PEO")
get.shared.areas(biovalue_dir="./3_8bio_extract", clade="subgen.Quercus",
                 stage_included=c("Ypresian", "Lutetian", "Bartonian", "Priabonian", "Rupelian", "Chattian"),
                 NorthHem_bio_dir="./1_rasterdata_NorthHem_asc_5min",
                 stage_mask=stage_mask_PEO, output_dir="./4_masked_bio_PEO")


