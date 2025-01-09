#####transfer modern coordinates to paleocoordinates#####
#https://gws.gplates.org/ 
library(chronosphere)


data_raw <- read.csv("../Fagaceae_Fossil_Records.csv")
colnames(data_raw)
genus_list <- c("Castanea", "Castanopsis", "Lithocarpus", "Chrysolepis", "Notholithocarpus", "Quercus")
data1 <- data_raw[data_raw$Genus%in%genus_list,]
data <- data1[!is.na(data1$Age_mean)&!is.na(data1$Longitude),]

data_cor <- data[,c("Longitude","Latitude")]
#age:(numeric)is the age in Ma at which the points will be reconstructed
age <- round(rowMeans(data[,c("Age_max","Age_min")]))
paleo_cor <- matrix(data=NA,nrow=nrow(data),ncol=2)

for (i in 1:nrow(data)){
  print(i)
  paleo_cor[i,] <- reconstruct(data_cor[i,], age[i])
}

paleo_cor <- as.data.frame(paleo_cor)
colnames(paleo_cor) <- c("paleoLong","paleoLat")
data<-cbind(data[,1:13],paleo_cor,age)
write.csv(data,file="FossilRecords_paleocor.csv")


