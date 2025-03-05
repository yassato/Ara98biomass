library(tidyverse)

###############################
# reshape SNPs for Wuest et al. (2022)

setwd("Atha/Wuest2022PLoSBiolData/")

geno = read.csv("sub_snps.csv.gz")
position = read.csv("positions.csv.gz")
colnames(geno)
geno = geno[,-1]
AF = apply(geno,1,mean)
hist(AF[which(AF>0.05&AF<0.95)],xlim=c(0,1))
geno = geno[which(AF>0.05&AF<0.95),]
position = position[which(AF>0.05&AF<0.95),]
MAF = AF[which(AF>0.05&AF<0.95)]
MAF[MAF>0.5] = 1 - MAF[MAF>0.5]

accs = read.csv("accs.csv")

accs$X0 = gsub("b","",accs$X0)
accs$X0 = gsub("'","",accs$X0)
accs$X0 = as.numeric(accs$X0)

gwasID = read.csv("gwasIDlist.csv")
colnames(geno) = gwasID$GWASid
saveRDS(geno,file="./SEW2022_sub_snpsMAF5.rds",compress=TRUE)

each = which((position$X0[2:length(position$X0)] - position$X0[1:(length(position$X0)-1)])<0)
chr = rep(c(1:5),times=c(each[1],each[2]-each[1],each[3]-each[2],each[4]-each[3],length(position$X0)-each[4]))
res = data.frame(chr,position$X0,MAF)
colnames(res) = c("chr","pos","maf")
saveRDS(res,file="./SEW_sub_posMAF5.rds",compress=TRUE)

