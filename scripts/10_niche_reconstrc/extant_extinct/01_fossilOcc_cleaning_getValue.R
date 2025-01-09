#2023-07-27 in R by LIUSY
library(chronosphere)
library(sp)
library(raster)

occ_df = read.csv("./20230727_66fossilspecies_occ.csv", header=T)

################################################################################
# (1) convert the present coordinates to paleo-coordinates
occ_df2 <- occ_df[!is.na(occ_df$Age_mean)&!is.na(occ_df$Longitude),]
data_cor <- occ_df2[,c("Longitude","Latitude")]

#age:(numeric)is the age in Ma at which the points will be reconstructed
age <- round(rowMeans(occ_df2[,c("Age_max","Age_min")]))
paleo_cor <- matrix(data=NA,nrow=nrow(occ_df2),ncol=2)

for (i in 1:nrow(occ_df2)){
  print(i)
  paleo_cor[i,] <- reconstruct(data_cor[i,], age[i])
}

paleo_cor <- as.data.frame(paleo_cor)
colnames(paleo_cor) <- c("paleoLong","paleoLat")
occ_df2 <- cbind(occ_df2[,1:10], paleo_cor, age)



################################################################################
# (2) find out the corresponding paleo-climate layer for each occurrence
layers_df = read.csv("./teyd_simulations.csv", header=T)

#assign the climate layer of closest geological age to each occurrence 
colnames(occ_df2)
occ_df2$layer = vector(length=nrow(occ_df2))

for (i in 1:nrow(occ_df2)){
  print(i)
  occ_df2$layer[i] = layers_df$Simulation_Name[which.min(abs(layers_df$Time-occ_df2$Age_mean[i]))]
}
#check
table(occ_df2$layer)



################################################################################
# (3) extract eight representative paleo-climate values for each occurrence
occ_df2$cellNumber = vector(length=nrow(occ_df2))
occ_df2$CMM = vector(length=nrow(occ_df2))
occ_df2$DryMon = vector(length=nrow(occ_df2))
occ_df2$MAP = vector(length=nrow(occ_df2))
occ_df2$MAT = vector(length=nrow(occ_df2))
occ_df2$WetDryMon = vector(length=nrow(occ_df2))
occ_df2$WetMon = vector(length=nrow(occ_df2))
occ_df2$WMM = vector(length=nrow(occ_df2))
occ_df2$WMMCMM = vector(length=nrow(occ_df2))

names <- c("CMM","DryMon","MAP","MAT","WetDryMon","WetMon","WMM","WMMCMM")

paleoclimate = unique(occ_df2$layer)
for (i in 1:length(paleoclimate)){
  print(paleoclimate[i])
  paleobio_dir = paste0("./teyd_asc_5min/",paleoclimate[i])
  files <- list.files(path =paleobio_dir, pattern = 'asc$|tif$',full.names=T)
  teu <- stack(files)
  print(teu)
  teu <- subset(teu, c(1,3,5,6,8,10,12,13)) #choose 8 monthly bio,for paleo-bio
  
  index = occ_df2$layer==paleoclimate[i]
  occ_df2$cellNumber[index] <- cellFromXY(teu[[1]], occ_df2[index,c('paleoLong','paleoLat')])
  
  sp_occ <- occ_df2[index,c('paleoLong','paleoLat')]
  coordinates(sp_occ) <- c('paleoLong','paleoLat')
  
  # extract bio climate and write it to your direction
  #If df=TRUE, the results are returned as a 'dataframe'##?CRS is NA
  occ_df2[index,names] = extract(teu, sp_occ, method='simple', buffer=1000, fun=mean, na.rm=TRUE, df=TRUE)[,2:9]
}

save(occ_df2,file="outputs/occ_df2_climates.Rdata")
num_na = occ_df2$WMMCMM[is.na(occ_df2$WMMCMM)]




extract_buffer <- function(df=occ_df4, buffer=1000000){
  for (i in 1:nrow(df)){
    print(i)
    if (is.na(df$WMMCMM[i])){
      paleobio_dir = paste0("./teyd_asc_5min/",df$layer[i])
      files <- list.files(path =paleobio_dir, pattern = 'asc$|tif$',full.names=T)
      teu <- stack(files)
      print(teu)
      teu <- subset(teu, c(1,3,5,6,8,10,12,13)) #choose 8 monthly bio,for paleo-bio
      
      sp_occ <- df[i,c('paleoLong','paleoLat')]
      coordinates(sp_occ) <- c('paleoLong','paleoLat')
      
      df[i,names] = extract(teu, sp_occ, method='simple', buffer=buffer, fun=mean, na.rm=TRUE, df=TRUE)[,2:9]
      print(df[i,names])
    }
  }
  return(df)
}

#for 22 NA-value, we extend the buffer of extracting value.
occ_df3 <- occ_df2
occ_df3 = extract_buffer(df=occ_df3, buffer=100000)
num_na2 = occ_df3$WMMCMM[is.na(occ_df3$WMMCMM)]
save(occ_df3,file="outputs/occ_df3_climates.Rdata")


#for 10 NA-value, we extend the buffer of extracting value.
occ_df4 <- occ_df3
occ_df4 = extract_buffer(df=occ_df4, buffer=1000000)
num_na3 = occ_df4$WMMCMM[is.na(occ_df4$WMMCMM)]
save(occ_df4,file="outputs/occ_df4_climates.Rdata")

# all occurrences with paleo-climate now


################################################################################
# (4) check occurrences within the same grid, and randomly retain a single occurrence
# i.e., 
load("outputs/occ_df4_climates.Rdata")

fossil_time_cell = paste(occ_df4$NicheRecon_tip,occ_df4$layer,occ_df4$cellNumber,sep="_")
dup = duplicated(fossil_time_cell)
table(dup)
occ_df4$dup.fossil_time_cell = dup
save(occ_df4,file="outputs/occ_df4_climates.Rdata")


write.csv(occ_df4, file="outputs/Fossil_paleocor_climates.csv",row.names = F,quote = F)






