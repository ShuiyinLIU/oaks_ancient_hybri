# This script is primarily modified from https://github.com/SheaML/ExDFOIL

library(ape)
library(phytools)
library(RColorBrewer)
library(ggtree)
library(ggimage)

setwd("E:/liushuiyin/Dfoil/HYB_98")
#setwd("E:/liushuiyin/Dfoil/RNA_2821")
###########################################################
# 1\ Getting DFOIL results into R
#now ready to read the table containing test results and associated sample information
#(all.summary.alt_appended.txt) into R.
MyDFOIL<-read.table("all.summary.alt_appended.txt",header=FALSE)
MyDFOIL<-as.data.frame(MyDFOIL)
# colnames(MyDFOIL)<-c("sp1","sp2","sp3","sp4","count1","count2","count3","count4",
#                      "introg","DFO_stat","DFO_p","DIL_stat","DIL_p","DFI_stat","DFI_p",
#                      "DOL_stat","DOL_p","name1","genus1","section1",
#                      "name2","genus2","section2","name3","genus3","section3",
#                      "name4","genus4","section4")
colnames(MyDFOIL)<-c("sp1","sp2","sp3","sp4","count1","count2","count3","count4",
                     "introg","DFO_stat","DFO_p","DIL_stat","DIL_p","DFI_stat","DFI_p",
                     "DOL_stat","DOL_p",
                     "name1","genus1","subgenus1","section1","region1",
                     "name2","genus2","subgenus2","section2","region2",
                     "name3","genus3","subgenus3","section3","region3",
                     "name4","genus4","subgenus4","section4","region4")
table(MyDFOIL$introg) #LSY
#RNA2821
#123    124     13     14     23     24     31     32     41     42   none 
#4161   3275    115     48    120     71     77     90     28     73 216940
#Proportion of positive tests 
sum(table(MyDFOIL$introg))-table(MyDFOIL$introg)["none"]
sum(table(MyDFOIL$introg))-table(MyDFOIL$introg)["none"]-table(MyDFOIL$introg)["na"]
#8058
(sum(table(MyDFOIL$introg))-table(MyDFOIL$introg)["none"])/sum(table(MyDFOIL$introg))
(sum(table(MyDFOIL$introg))-table(MyDFOIL$introg)["none"]-table(MyDFOIL$introg)["na"])/sum(table(MyDFOIL$introg))
#0.03581365 for RNA2821
#Proportion of “ancestral” signatures
(table(MyDFOIL$introg)["123"]+table(MyDFOIL$introg)["124"])/sum(table(MyDFOIL$introg))
#0.03304918 for RNA2821
#Proportion of “intergroup” signatures
(sum(table(MyDFOIL$introg))-table(MyDFOIL$introg)["none"]-table(MyDFOIL$introg)["123"]-table(MyDFOIL$introg)["124"])/sum(table(MyDFOIL$introg))
#0.002764469 for RNA2821

###########################################################
# 2\ Classifying introgression
#helpful to bin introgression signatures into categories or "classes". 
#introg_classifier.R that will translate the "raw" DFOIL introgression result
#(i.e., the format returned by dfoil_analyze.py; e.g., "123" or "31") into an 
#introgression class, which might be based on species, subspecies, locality, or
#any other column of the sample information data frame. 

#For reusability, the function requires that you specify the column number containing
#raw DFOIL introgression results, as well as the column numbers corresponding to the
#categories you wish to bin by for each taxon P1 through P4. I use apply() to execute
#this function on each row of the dataframe, with column numbers passed to apply after
#the function name.
source("../introg_classifier.R")
#firstly summarize results between genera
#9: column "introg"; 19: column "genus1"; 24: column "genus2"; 29: column "genus3"; 34: column "genus4"
MyDFOIL$introg_class <- unlist(apply(MyDFOIL,1,introg_classifier,9,19,24,29,34))

###########################################################
# 3\ Summary tables
#return the number of test results in each introgression "class"
table(MyDFOIL$introg_class) 

#The percent of each introgression “class" in total positive introgression tests
table(MyDFOIL[(MyDFOIL$introg_class != "none")&(MyDFOIL$introg_class != "na"),]$introg_class)
prop.table(table(MyDFOIL[(MyDFOIL$introg_class != "none")&(MyDFOIL$introg_class != "na"),]$introg_class))
#class "na" will be included into class “none“
write.csv(table(MyDFOIL$introg_class),"introg_class_genus_df.csv")

###########################################################
# 4\ Summary and visualization of results across phylogeny
#To summarize DFOIL results across phylogenetic space, we will first need read in a 
#rooted phylogenetic tree containing all of the individuals. This tree can optionally
#be pruned to remove individuals not involved in any tests. 
#MyTree <- read.tree("mcctree.conRT89_2821k_top20_10000resample_74dfoil.newick")
MyTree <- read.tree("mcctree.calib0a9f_conMO431_89k_top20_10000resample_150dfoil.newick")


#use phytools getMRCA() function to return the node number of the most recent 
#common ancestor of taxa P1 and P2 ("mrca12") and taxa P3 and P4 ("mrca34").
MyDFOIL$mrca12<-apply(MyDFOIL,1,function(x) getMRCA(MyTree,c(x[1],x[2])))
MyDFOIL$mrca34<-apply(MyDFOIL,1,function(x) getMRCA(MyTree,c(x[3],x[4])))
MyDFOIL$mrca12 -> mrca12
MyDFOIL$mrca34 -> mrca34
save(mrca12, file = "mrca12.RData")
save(mrca34, file = "mrca34.RData")
#load("mrca12.RData")
#load("mrca34.RData")
#MyDFOIL$mrca12 <- mrca12
#MyDFOIL$mrca34 <- mrca34

#Now that the node numbers are associated with each test, we can summarize test 
#results by the ancestral node(s) involved again using aggregate().
aggregate(as.factor(mrca12) ~ introg_class,data=MyDFOIL,FUN=table) -> mrca12_fig
aggregate(as.factor(mrca34) ~ introg_class,data=MyDFOIL,FUN=table) -> mrca34_fig

save(mrca12_fig, file = "mrca12_fig_genus.RData")
save(mrca34_fig, file = "mrca34_fig_genus.RData")
#load("mrca12_fig_genus.RData")
#load("mrca34_fig_genus.RData")


#Next, we reformat the output for ggtree:
#For the MRCA of Taxa P1 and P2
as.matrix(mrca12_fig[,-1]) -> mrca12_fig_num 
rownames(mrca12_fig_num) <- mrca12_fig[,1]
prop.table(mrca12_fig_num) -> mrca12_fig_prop
t(mrca12_fig_prop) -> mrca12_fig_prop_t
as.data.frame(mrca12_fig_prop_t) -> mrca12_fig_prop_df
mrca12_fig_prop_df$node <- rownames(mrca12_fig_prop_df)

#For the MRCA of Taxa P3 and P4
as.matrix(mrca34_fig[,-1]) -> mrca34_fig_num 
rownames(mrca34_fig_num) <- mrca34_fig[,1]
prop.table(mrca34_fig_num) -> mrca34_fig_prop
t(mrca34_fig_prop) -> mrca34_fig_prop_t
as.data.frame(mrca34_fig_prop_t) -> mrca34_fig_prop_df
mrca34_fig_prop_df$node <- rownames(mrca34_fig_prop_df)

# set color for each introgression class
#
#col_introgenus=data.frame(introg_class=colnames(mrca12_fig_prop_df)[-38],
#                          introg_ID=c(rep("Castanea -> X",3),
#                                      rep("Castanea / Castanea <-> X",4),
#                                      rep("Castanopsis -> X",2),
#                                      rep("Castanopsis / Castanea <-> X",4),
#                                     rep("Castanopsis / Castanopsis <-> X",3),
#                                      rep("Lithocarpus -> X",5),
#                                      rep("Lithocarpus / Lithocarpus <-> X",5),
#                                      "none",
#                                      rep("Notholithocarpus -> X",2),
#                                      rep("Quercus -> X",4),
#                                      rep("Quercus / Quercus <-> X",4)),
#                          cols=c(rep("orange",3),
#                                 rep("yellow",4),
#                                 rep("purple",2),
#                                 rep("brown",4),
#                                 rep("red",3),
#                                 rep("lightgreen",5),
#                                 rep("green",5),
#                                 "gray",
#                                 rep("pink",2),
#                                 rep("blue",4),
#                                 rep("cyan",4)))
##ONLY introg class related Quercus is given a certain color；this for HYB-98
others_c="black"
col_introgenus=data.frame(introg_class=colnames(mrca12_fig_prop_df)[-36],
                          introg_ID=c("others","Castanea -> Quercus",
                                      rep("others",2),"Castanea / Castanea <-> Quercus",
                                      "others","Castanopsis -> Quercus",
                                      "others","Castanopsis / Castanea <-> Quercus",
                                      rep("others",3),
                                      "Castanopsis / Castanopsis <-> Quercus",
                                      rep("others",3),
                                      "Chrysolepis / Chrysolepis <-> Quercus",
                                      rep("others",3),
                                      "Lithocarpus -> Quercus",
                                      rep("others",4),
                                      "Lithocarpus / Lithocarpus <-> Quercus",
                                      "none",
                                      "Quercus -> Castanea",
                                      "Quercus -> Castanopsis",
                                      "Quercus -> Lithocarpus",
                                      "Quercus -> Quercus",
                                      "Quercus / Quercus <-> Castanea",
                                      "Quercus / Quercus <-> Castanopsis",
                                      "Quercus / Quercus <-> Lithocarpus",
                                      "Quercus / Quercus <-> Quercus"),
                          cols=c(others_c,"Khaki",rep(others_c,2),"yellow",
                                 others_c,"orange",others_c,"DarkGoldenrod",
                                 rep(others_c,3),"pink",rep(others_c,3),"Coral",
                                 rep(others_c,3),"red",rep(others_c,4),"brown",
                                 "gray","Orchid","purple","lightgreen","SpringGreen",
                                 "SeaGreen","cyan","SteelBlue","blue"))
##ONLY introg class related Quercus is given a certain color；this for RNA-2821
#others_c="black"
#col_introgenus=data.frame(introg_class=colnames(mrca12_fig_prop_df)[-38],
#                          introg_ID=c(rep("others",2),"Castanea -> Quercus",
#                                      rep("others",3),"Castanea / Castanea <-> Quercus",
#                                      "others","Castanopsis -> Quercus",
#                                      rep("others",3),"Castanopsis / Castanea <-> Quercus",
#                                      rep("others",2),
#                                      "Castanopsis / Castanopsis <-> Quercus",
#                                      rep("others",4),
#                                      "Lithocarpus -> Quercus",
#                                      rep("others",4),
#                                      "Lithocarpus / Lithocarpus <-> Quercus",
#                                      "none",
#                                      rep("others",2),
#                                      "Quercus -> Castanea",
#                                      "Quercus -> Castanopsis",
#                                      "Quercus -> Lithocarpus",
#                                      "Quercus -> Quercus",
#                                      "Quercus / Quercus <-> Castanea",
#                                      "Quercus / Quercus <-> Castanopsis",
#                                      "Quercus / Quercus <-> Lithocarpus",
#                                      "Quercus / Quercus <-> Quercus"),
#                          cols=c(rep(others_c,2),"Khaki",rep(others_c,3),"yellow",
#                                 others_c,"orange",rep(others_c,3),"DarkGoldenrod",
#                                 rep(others_c,2),"pink",rep(others_c,4),
#                                 "red",rep(others_c,4),"brown","gray",rep(others_c,2),
#                                 "Orchid","purple","lightgreen","SpringGreen",
#                                 "SeaGreen","cyan","SteelBlue","blue"))

#plot for the MRCA of Taxa P1 and P2
pies12<-nodepie(mrca12_fig_prop_df,cols=1:35,color = col_introgenus$cols)

#plot for the MRCA of Taxa P3 and P4
pies34<-nodepie(mrca34_fig_prop_df,cols=1:35,color = col_introgenus$cols)

#Finally, for plotting with ggtree:
MyTree_GG <- ggtree(MyTree)+theme_tree2() 
pietest<-inset(MyTree_GG + geom_tiplab(size=1.8)+ ggplot2::xlim(0,75)+ ggplot2::ylim(1,150),pies34, height=0.06, width=0.06, hjust=1.4,vjust=-1.3)
pietest2<-inset(pietest,pies12, height=0.06, width=0.06, hjust=1.4,vjust=1.9)

pdf("intro_class_pie_genuslevel_timescale.pdf",width = 9,height = 11)
print(pietest2)
dev.off()
#add legend
pdf("intro_class_pie_genuslevel_legend.pdf",width = 9,height = 10)
plot(0,0,type="n",ann = F, bty = "n",xaxt = "n", yaxt = "n")
col_introgenus_uniqe <- col_introgenus[!duplicated(col_introgenus$introg_ID),]
legend(-1,1,legend=col_introgenus_uniqe$introg_ID,col=col_introgenus_uniqe$cols,pch=16,cex=1)
dev.off()



###############################################################################
###############################################################################
#at section levels
colnames(MyDFOIL)
MyDFOIL$introg_class_sect <- unlist(apply(MyDFOIL,1,introg_classifier,9,21,26,31,36))
table(MyDFOIL$introg_class_sect) 

table(MyDFOIL[MyDFOIL$introg_class_sect != "none",]$introg_class_sect)
prop.table(table(MyDFOIL[MyDFOIL$introg_class_sect != "none",]$introg_class_sect))

write.csv(table(MyDFOIL$introg_class_sect),"introg_class_sect_df.csv")
#extract and summary the total tests number of each introg_class
introg_class_sect_df <- read.csv("introg_class_sect_df.csv",header=T,stringsAsFactors = F)
intro_tests <- vector(length=nrow(introg_class_sect_df))
for (i in 1:nrow(introg_class_sect_df)){
  print(i)
  if (introg_class_sect_df$tests_class[i]=="and"){
    intro_tests[i] <- nrow(MyDFOIL[(MyDFOIL$section1==introg_class_sect_df$P1[i]&MyDFOIL$section2==introg_class_sect_df$P2[i])&
      (MyDFOIL$section3==introg_class_sect_df$P34[i]|MyDFOIL$section4==introg_class_sect_df$P34[i]),])
  }else
    intro_tests[i] <- nrow(MyDFOIL[(MyDFOIL$section1==introg_class_sect_df$P1[i]|MyDFOIL$section2==introg_class_sect_df$P1[i]|
                                    MyDFOIL$section3==introg_class_sect_df$P1[i]|MyDFOIL$section4==introg_class_sect_df$P1[i])&
                                     (MyDFOIL$section1==introg_class_sect_df$P34[i]|MyDFOIL$section2==introg_class_sect_df$P34[i]|
                                        MyDFOIL$section3==introg_class_sect_df$P34[i]|MyDFOIL$section4==introg_class_sect_df$P34[i]),])
}
introg_class_sect_df$intro_tests <- intro_tests
write.csv(introg_class_sect_df,"introg_class_sect_df.csv")

#if ("CY"%in%sect_p1234[1]){print("in")}else{print("no")}
#library(tidyverse)
#grepl("CY",sect_p1234[1:10])&grepl("CY",sect_p1234[1:10])

###############################################################################
###############################################################################
#at subgenus levels
source("../introg_classifier.R")
colnames(MyDFOIL)
MyDFOIL$introg_class_subgenus <- unlist(apply(MyDFOIL,1,introg_classifier,9,20,25,30,35))
table(MyDFOIL$introg_class_subgenus) 

table(MyDFOIL[MyDFOIL$introg_class_subgenus != "none",]$introg_class_subgenus)
prop.table(table(MyDFOIL[MyDFOIL$introg_class_subgenus != "none",]$introg_class_subgenus))

write.csv(table(MyDFOIL$introg_class_subgenus),"introg_class_subgenus_df.csv")



###########################################################
###########################################################
##at regional levels (old world [OW] and new world [NW])
source("../introg_classifier.R")
colnames(MyDFOIL)
MyDFOIL$introg_class_region <- unlist(apply(MyDFOIL,1,introg_classifier,9,22,27,32,37))
table(MyDFOIL$introg_class_region) 

table(MyDFOIL[MyDFOIL$introg_class_region != "none",]$introg_class_region)
prop.table(table(MyDFOIL[MyDFOIL$introg_class_region != "none",]$introg_class_region))

write.csv(table(MyDFOIL$introg_class_region),"introg_class_region_df.csv")






