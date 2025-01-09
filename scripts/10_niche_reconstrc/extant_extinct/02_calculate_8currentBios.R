#created at 2023-07-28 15:17 in R by LIUSY

library(raster)


files <- list.files(path ="../wc2.1_30s_tavg/", pattern = 'asc$|tif$',full.names=T)
s <- stack(files)
print(s)




tempMon_max = calc(s, max)
tempMon_min = calc(s, min)
tempMon_max_min_diff = overlay(tempMon_max, tempMon_min, fun=function(x,y){return(x-y)})

writeRaster(tempMon_max, "wc2.1_30s_tavg.wmm.tif", overwrite=TRUE)
writeRaster(tempMon_min, "wc2.1_30s_tavg.cmm.tif", overwrite=TRUE)
writeRaster(tempMon_max_min_diff, "wc2.1_30s_tavg.wmmcmm.tif", overwrite=TRUE)


prepMon_max = raster("/nfs/liushuiyin/20211115_oakhyb_add/niche_reconstrc/EnvLayer/wc2.0_bio_30s_13.asc")
prepMon_min = raster("/nfs/liushuiyin/20211115_oakhyb_add/niche_reconstrc/EnvLayer/wc2.0_bio_30s_14.asc")  
prepMon_max_min_diff = overlay(prepMon_max, prepMon_min, fun=function(x,y){return(x-y)})
writeRaster(prepMon_max_min_diff, "wc2.0_bio_30s.wetdrymon.tif", overwrite=TRUE)





  