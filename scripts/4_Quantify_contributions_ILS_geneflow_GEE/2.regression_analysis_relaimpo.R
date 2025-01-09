#####
#compile data
library(ape)
library(phytools)

##### RNA2821
setwd("E:/liushuiyin/CS_FGdetection/RT89_2821")
astral_sp = read.tree("astral_RT89_MOot.V1_2821.tre.rr")
tr_CF_raw = read.tree("concord.cf.tree")
tr_ERR_raw = read.tree("RAxML_bipartitions.ERR.rr")
tr_mlfixed_raw = read.tree("RAxML_bestTree.astralconstrained.tre.rr")
tr_RI_raw = read.tree("unbalanced_triples_perc_reticulation_index.tre")

trees <- c(astral_sp, tr_CF_raw, tr_ERR_raw, tr_mlfixed_raw, tr_RI_raw) 
trees <- .compressTipLabel(trees)
tr_CF <- trees[[2]]
tr_ERR <- trees[[3]]
tr_mlfixed <- trees[[4]]
tr_RI <- trees[[5]]

all.equal(astral_sp$tip.label, tr_CF$tip.label)

#######################################################################
##### HYB98
setwd("E:/liushuiyin/CS_FGdetection/RT431_98")
astral_sp = read.tree("astral_RT_ot8ot7ot2R2.V3sra.tm_431.V1_98.tre.rr.rr")
tr_CF_raw = read.tree("concord.cf.tree.rr")
tr_ERR_raw = read.tree("RAxML_bipartitions.ERR")
tr_mlfixed_raw = read.tree("RAxML_bestTree.astralconstrained.tre.rr")
tr_RI_raw = read.tree("unbalanced_triples_perc_reticulation_index.tre.rr")

trees <- c(astral_sp, tr_CF_raw, tr_ERR_raw, tr_mlfixed_raw, tr_RI_raw) 
trees <- .compressTipLabel(trees)
tr_CF <- trees[[2]]
tr_ERR <- trees[[3]]
tr_mlfixed <- trees[[4]]
tr_RI <- trees[[5]]

all.equal(astral_sp$tip.label, tr_CF$tip.label)



#######################################################################
#######################################################################

#check the node number of all trees if they are same
df=data.frame(tr_CF=tr_CF$tip.label,tr_ERR$tip.label,astral_sp$tip.label,tr_mlfixed$tip.label,tr_RI=tr_RI$tip.label)
pdf("tr_CF.nodenumber.pdf",width=8,height=10)
plot(tr_CF,cex=0.5)
nodelabels(cex=0.7)
dev.off()
pdf("tr_ERR.nodenumber.pdf",width=8,height=10)
plot(tr_ERR,cex=0.5)
nodelabels(cex=0.7)
dev.off()
pdf("astral_sp.nodenumber.pdf",width=8,height=10)
plot(astral_sp,cex=0.5)
nodelabels(cex=0.7)
dev.off()
pdf("tr_mlfixed.nodenumber.pdf",width=8,height=10)
plot(tr_mlfixed,cex=0.5)
nodelabels(cex=0.7)
dev.off()
pdf("tr_RI.nodenumber.pdf",width=8,height=10)
plot(tr_RI,cex=0.5)
nodelabels(cex=0.7)
dev.off()


#create a dataframe
nodeID <- (length(astral_sp$tip.label)+1):(length(astral_sp$tip.label)+astral_sp$Nnode)
Y_BSsum = vector(length=length(nodeID))
for (i in 1:length(tr_CF$node.label)){Y_BSsum[i] <- as.numeric(unlist(strsplit(as.character(tr_CF$node.label[i]), split = "/")))[2]}
X3_Hyb <- tr_RI$node.label
X2_Err <- tr_ERR$node.label
#theta = mutation units/coalescent units
X1_ILS <- vector(length=length(nodeID))
for (i in 1:length(nodeID)){
  if (nodeID[i]%in%astral_sp$edge[,2]){
    X1_ILS[i] <- tr_mlfixed$edge.length[which(tr_mlfixed$edge[,2]==nodeID[i])]/
      astral_sp$edge.length[which(astral_sp$edge[,2]==nodeID[i])]
  }else{
    X1_ILS[i] <- NA
    }
  }

RC_df <- data.frame(nodeID=nodeID,Y_BSsum=Y_BSsum,X3_Hyb=X3_Hyb,X2_Err=X2_Err,X1_ILS=X1_ILS)
write.csv(RC_df,"relative_contribution_ILS_Err_Intro.csv")



################################################################################
library(relaimpo)
#This R script conducts relative importance decomposition to assign shares of importance of each variable and plot the results.

x_raw=read.csv('relative_contribution_ILS_Err_Intro.csv')
x=x_raw[(x_raw$X.1=="IN")&(x_raw$X1_ILS!=Inf),]


#without log transformations
data=x[,c('Y_BSsum','X1_ILS','X2_Err','X3_Hyb')]
bootimpo.result <- boot.relimp(data, b = 100,
                    type = c("lmg", "last", "first", "pratt"),
                    rank = TRUE, diff = TRUE, rela = TRUE)
booteval.relimp(bootimpo.result)

##Adding interaction terms in the regression model 
linmod <- lm(Y_BSsum ~ log(X1_ILS)+log(X2_Err)+log(X3_Hyb)+log(X1_ILS * X2_Err) + log(X1_ILS * X3_Hyb) + log(X2_Err * X3_Hyb)+log(X1_ILS * X3_Hyb*X2_Err), data = data)
summary(linmod)

## Plot
pdf("relative_importances_1Fagaceae2oaks.pdf",width=10,height=8)
plot(booteval.relimp(bootimpo.result))

##### only genus Quercus
x=x_raw[(x_raw$X.1=="IN")&(x_raw$X1_ILS!=Inf)&(x_raw$X.2=="oak"),]
data=x[,c('Y_BSsum','X1_ILS','X2_Err','X3_Hyb')]
bootimpo.result <- boot.relimp(data, b = 100,
                               type = c("lmg", "last", "first", "pratt"),
                               rank = TRUE, diff = TRUE, rela = TRUE)
booteval.relimp(bootimpo.result)
##Adding interaction terms in the regression model 
linmod <- lm(Y_BSsum ~ log(X1_ILS)+log(X2_Err)+log(X3_Hyb)+log(X1_ILS * X2_Err) + log(X1_ILS * X3_Hyb) + log(X2_Err * X3_Hyb)+log(X1_ILS * X3_Hyb*X2_Err), data = data)
summary(linmod)

## Plot
plot(booteval.relimp(bootimpo.result))

dev.off()



#####################################################
## apply log transformations to the regressors
x=x_raw[(x_raw$X.1=="IN")&(x_raw$X1_ILS!=Inf),]
data=x[,c('Y_BSsum','X1_ILS','X2_Err','X3_Hyb')]
data$X1_ILS <- log(data$X1_ILS)
data$X2_Err <- log(data$X2_Err)
data$X3_Hyb <- log(data$X3_Hyb)
data2 <- data[(data$X1_ILS!=Inf)&(data$X2_Err!=Inf)&(data$X3_Hyb!=Inf)&
                (data$X1_ILS!=-Inf)&(data$X2_Err!=-Inf)&(data$X3_Hyb!=-Inf),]
bootimpo.result <- boot.relimp(data2, b = 100,
                               type = c("lmg", "last", "first", "pratt"),
                               rank = TRUE, diff = TRUE, rela = TRUE)
booteval.relimp(bootimpo.result)
##Adding interaction terms in the regression model 
linmod <- lm(Y_BSsum ~ log(X1_ILS)+log(X2_Err)+log(X3_Hyb)+log(X1_ILS * X2_Err) + log(X1_ILS * X3_Hyb) + log(X2_Err * X3_Hyb)+log(X1_ILS * X3_Hyb*X2_Err), data = data2)
summary(linmod)

## Plot
pdf("relative_importances_1Fagaceae2oaks_logtranf.pdf",width=10,height=8)
plot(booteval.relimp(bootimpo.result))

##### only genus Quercus
x=x_raw[(x_raw$X.1=="IN")&(x_raw$X1_ILS!=Inf)&(x_raw$X.2=="oak"),]
data=x[,c('Y_BSsum','X1_ILS','X2_Err','X3_Hyb')]
data$X1_ILS <- log(data$X1_ILS)
data$X2_Err <- log(data$X2_Err)
data$X3_Hyb <- log(data$X3_Hyb)
data2 <- data[(data$X1_ILS!=Inf)&(data$X2_Err!=Inf)&(data$X3_Hyb!=Inf)&
                (data$X1_ILS!=-Inf)&(data$X2_Err!=-Inf)&(data$X3_Hyb!=-Inf),]
bootimpo.result <- boot.relimp(data2, b = 100,
                               type = c("lmg", "last", "first", "pratt"),
                               rank = TRUE, diff = TRUE, rela = TRUE)
booteval.relimp(bootimpo.result)
##Adding interaction terms in the regression model 
linmod <- lm(Y_BSsum ~ log(X1_ILS)+log(X2_Err)+log(X3_Hyb)+log(X1_ILS * X2_Err) + log(X1_ILS * X3_Hyb) + log(X2_Err * X3_Hyb)+log(X1_ILS * X3_Hyb*X2_Err), data = data2)
summary(linmod)

## Plot
plot(booteval.relimp(bootimpo.result))

dev.off()

