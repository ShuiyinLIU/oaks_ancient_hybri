
library(raster)


crop.paleobio <- function(paleobio_dir, stage="Thanetian", output_dir){
  files <- list.files(path =paleobio_dir, pattern = 'asc$',full.names=T)
  teu <- stack(files)
  print(teu)
  teu <- subset(teu, c(1,3,5,6,8,10,12,13)) #choose 8 monthly bio,for paleo-bio
  print(teu)
  # extract extent
  bios <- crop(teu, extent(-180,180,-20,90)) #mainly North Hemisphere
  print(names(bios))
  teu_name<-as.list(paste(stage, c("NorthHem_cmm","NorthHem_drymon","NorthHem_map",
                 "NorthHem_mat","NorthHem_wetdrymon","NorthHem_wetmon",
                 "NorthHem_wmm","NorthHem_wmmcmm"), sep="_")) #rename,for paleo-bio
  # write
  for (i in 1:8) {
    writeRaster(bios[[i]], filename=paste(output_dir, teu_name[[i]], sep="/"), format="ascii", overwrite=TRUE)
  }
}

crop.paleobio(paleobio_dir="H:/teyd_mask_asc_5min/teydn1", stage="Danian",
              output_dir="./1_rasterdata_NorthHem_asc_5min")
crop.paleobio(paleobio_dir="H:/teyd_mask_asc_5min/teydm1", stage="Thanetian",
              output_dir="./1_rasterdata_NorthHem_asc_5min")
crop.paleobio(paleobio_dir="H:/teyd_mask_asc_5min/teydl1", stage="Ypresian",
              output_dir="./1_rasterdata_NorthHem_asc_5min")
crop.paleobio(paleobio_dir="H:/teyd_mask_asc_5min/teydj1", stage="Lutetian",
              output_dir="./1_rasterdata_NorthHem_asc_5min")
crop.paleobio(paleobio_dir="H:/teyd_mask_asc_5min/teydi1", stage="Bartonian",
              output_dir="./1_rasterdata_NorthHem_asc_5min")
crop.paleobio(paleobio_dir="H:/teyd_mask_asc_5min/teydh1", stage="Priabonian",
              output_dir="./1_rasterdata_NorthHem_asc_5min")
crop.paleobio(paleobio_dir="H:/teyd_mask_asc_5min/teydg1", stage="Rupelian",
              output_dir="./1_rasterdata_NorthHem_asc_5min")
crop.paleobio(paleobio_dir="H:/teyd_mask_asc_5min/teydf1", stage="Chattian",
              output_dir="./1_rasterdata_NorthHem_asc_5min")


# check
checks <- list.files(path = "./1_rasterdata_NorthHem_asc_5min", pattern = 'asc$',full.names=T)
check <- stack(checks)
check
plot(check,4)


