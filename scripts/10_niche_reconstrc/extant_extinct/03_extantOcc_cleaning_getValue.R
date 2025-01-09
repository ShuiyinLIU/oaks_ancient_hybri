#created at 2023-07-28 11:10 in R by LIUSY
library(raster)

geo_df_raw <- read.csv("D:/liushuiyin/20211110_hybridization/20230511投稿三稿/20230605RevisedV1/Newanalyses_with_fossils/03.nicheRecon/extant_only/01.geodata_cleaned_421sel.csv", header=T)
geo_name_df <- read.csv("D:/liushuiyin/20211110_hybridization/20230511投稿三稿/20230605RevisedV1/Newanalyses_with_fossils/03.nicheRecon/extant_only/05.geoSname_hybSname.csv", header=F)


##############################################################################################
# (1) process the occurrences
#only these species with more than three occurrences (i.e., 409 extant species) are used for extracting climate value and niche recon.
dim(geo_name_df)
geo_df <- data.frame(species=vector(),longitude=vector(),latitude=vector(),renameSequencingData=vector())

for (i in 1:nrow(geo_name_df)){
  print(i)
  df_sub = geo_df_raw[geo_df_raw$species==geo_name_df$V1[i],]
  df_sub$renameSequencingData = rep(geo_name_df$V2[i],nrow(df_sub))
  geo_df = rbind(geo_df, df_sub)
}

dim(geo_df)

# add a column for "dup.species_cell" for keep only one occurrence for one species in each cell
rs <- raster(list.files(path ="C:/Users/liush/Downloads/wc2.1_30s_tavg/", pattern = 'asc$|tif$',full.names=T)[1])
#rs <- raster(list.files(path ="./wc2.1_30s_tavg/", pattern = 'asc$|tif$',full.names=T)[1])

geo_df$cellNumber <- cellFromXY(rs, geo_df[,c('longitude','latitude')])

species_cell = paste(geo_df$renameSequencingData,geo_df$cellNumber,sep="_")
dup = duplicated(species_cell)
table(dup)
geo_df$dup.species_cell = dup




##############################################################################################
# (2) extract 8 current variables for each occurrence

geo_df$CMM = vector(length=nrow(geo_df))
geo_df$DryMon = vector(length=nrow(geo_df))
geo_df$MAP = vector(length=nrow(geo_df))
geo_df$MAT = vector(length=nrow(geo_df))
geo_df$WetDryMon = vector(length=nrow(geo_df))
geo_df$WetMon = vector(length=nrow(geo_df))
geo_df$WMM = vector(length=nrow(geo_df))
geo_df$WMMCMM = vector(length=nrow(geo_df))

save(geo_df,file="geo_df.Rdata")

load("./geo_df.Rdata")





extract_buffer <- function(df=geo_df, buffer=1000){
  names <- c("CMM","DryMon","MAP","MAT","WetDryMon","WetMon","WMM","WMMCMM")
  files <- list.files(path ="./currentWorldclim", pattern = 'asc$|tif$',full.names=T)
  teu <- stack(files)
  print(teu)
  #adjust the position of each varible to match the positon in fossil climateValue fossil df.
  teu <- subset(teu, c(6,4,2,1,5,3,7,8)) #choose 8 monthly bio
  print(teu)
  
  sp_occ <- df[,c('longitude','latitude')]
  coordinates(sp_occ) <- c('longitude','latitude')
  
  df[,names] = extract(teu, sp_occ, method='simple', buffer=buffer, fun=mean, na.rm=TRUE, df=TRUE)[,2:9]

  return(df)
}



geo_df = extract_buffer(df=geo_df, buffer=1000)
num_na = geo_df$WMMCMM[is.na(geo_df$WMMCMM)]

save(geo_df,file="outputs/geo_df_climatesCurrent.Rdata")

print(length(num_na))




load("outputs/geo_df_climatesCurrent.Rdata")








