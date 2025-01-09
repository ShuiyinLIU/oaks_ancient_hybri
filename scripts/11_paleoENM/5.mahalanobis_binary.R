#####calculate MD#####
#https://pvanb.wordpress.com/2014/05/13/a-new-method-and-tool-exdet-to-evaluate-novelty-environmental-conditions/
library(raster)
#library(RColorBrewer)
library(colorRamps)
# Names of the reference (ref) and projection (pro) data



calculate_MD <- function(ref_dir="./4_masked_bio_PE", clade="Castanea", stage="Bartonian",
                         pro_dir="./1_rasterdata_NorthHem_asc_5min", output_dir="./5_MD_raster"){
  # Import the data in R using the read.asciigrid function of the sp package
  ##ref, cropped and masked bio layers
  ref_prefix <- paste(ref_dir,paste(clade,stage,"NorthHem",sep="_"),sep="/")
  ref <- paste(ref_prefix,c('cmm.asc','drymon.asc','map.asc','mat.asc','wetdrymon.asc','wetmon.asc','wmm.asc','wmmcmm.asc'),sep="_")
  refdat <- list()
  for(i in 1:8){
    refdat[[i]] <- read.asciigrid(fname=ref[i])@data
  }
  ##pro, cropped yet not masked bio layers
  pro_prefix <- paste(pro_dir,paste(stage,"NorthHem",sep="_"),sep="/")
  pro <- paste(pro_prefix,c('cmm.asc','drymon.asc','map.asc','mat.asc','wetdrymon.asc','wetmon.asc','wmm.asc','wmmcmm.asc'),sep="_")
  prodat <- list()
  for(j in 1:8){
    prodat[[j]] <- read.asciigrid(fname=pro[j])@data
  }
  refdat <- do.call(cbind, refdat)
  prodat <- do.call(cbind, prodat)
  print("Import ref and pro data finshed!")
  
  # Calculate the average and covariance matrix of the variables in the reference set
  ref.av  <- colMeans(refdat, na.rm=TRUE)
  ref.cov <- var(refdat, na.rm=TRUE)
  print(ref.av)
  
  # Calculate the mahalanobis distance of each raster 
  # cell to the environmental center of the reference 
  # set for both the reference and the projection data 
  # set and calculate the ratio between the two.
  mah.ref <- mahalanobis(x=refdat, center=ref.av, cov=ref.cov, tol=1e-20)
  mah.pro <- mahalanobis(x=prodat, center=ref.av, cov=ref.cov, tol=1e-20)
  ref.max <- max(mah.ref[is.finite(mah.ref)])
  nt1 <- as.data.frame(mah.ref/ref.max)
  pro.max <- max(mah.pro[is.finite(mah.pro)])
  nt2 <- as.data.frame(mah.pro/pro.max)
  
  #calculate log of nt, so as to magnify the small number, but you may not need to do this
  nt1<-log(nt1)
  nt2<-log(nt2)
  nt2[!is.na(nt2)]<-nt2[!is.na(nt2)]-min(na.omit(nt2))#seperate the range of nt1 and nt2
  range(na.omit(nt1))
  range(na.omit(nt2))
  
  # Create and plot the raster layer
  ###plot extract range
  # NT1 <- read.asciigrid(fname=pro[1])
  # NT1@data <- nt1
  # NT1rast <- raster(NT1)
  # plot(NT1rast, col=rev(blue2green2red(100)),main='ref')
  
  NT2 <- read.asciigrid(fname=pro[1])
  NT2@data <- nt2
  NT2rast <- raster(NT2)
  plot(NT2rast, col=rev(blue2green2red(100)),main='pro')
  
  #rescale to 0-1
  rescale <- function(x){(x-cellStats(x, "min"))/(cellStats(x, "max")-cellStats(x, "min"))}
  NT2rast <- rescale(NT2rast)
  #plot(NT2rast, col=rev(blue2green2red(100)),main=paste(clade,stage, sep=": "))
  
  #save as raster files
  crs(NT2rast) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
  writeRaster(NT2rast, paste(output_dir,paste(clade,stage, sep="_"),sep="/"), format="ascii", overwrite=T)
  writeRaster(NT2rast, paste(output_dir,paste(clade,stage, sep="_"),sep="/"), format="GTiff", overwrite=T)
}


clade_list=c("Castanea","Castanopsis","Lithocarpus","subgen.Cerris","subgen.Quercus")
stage_mask_PE=c("Danian","Thanetian","Ypresian", "Lutetian", "Bartonian", "Priabonian")
stage_mask_PEO=c("Danian","Thanetian","Ypresian", "Lutetian", "Bartonian", "Priabonian", "Rupelian", "Chattian")


######## PE
for (i in 1:length(clade_list)){
  print(clade_list[i])
  for (j in 1:length(stage_mask_PE)){
    print(stage_mask_PE[j])
    calculate_MD(ref_dir="./4_masked_bio_PE", clade=clade_list[i], stage=stage_mask_PE[j],
                 pro_dir="./1_rasterdata_NorthHem_asc_5min", output_dir="./5_MD_raster_PE")
  }
}



######## PEO
for (i in 1:length(clade_list)){
  print(clade_list[i])
  for (j in 1:length(stage_mask_PEO)){
    print(stage_mask_PEO[j])
    calculate_MD(ref_dir="./4_masked_bio_PEO", clade=clade_list[i], stage=stage_mask_PEO[j],
                 pro_dir="./1_rasterdata_NorthHem_asc_5min", output_dir="./5_MD_raster_PEO")
  }
}






