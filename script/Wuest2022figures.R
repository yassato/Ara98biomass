library(tidyverse)
library(patchwork)
library(rNeighborGWAS)
library(gaston)
source("coord.R")

############
# Figure 1
fA = function(x,b0,b1,b2,b12) { return((b12+b2)*(2*x-1)+b0+b1) }
fa = function(x,b0,b1,b2,b12) { return((b12-b2)*(2*x-1)+b0-b1) }

W_A = function(x,b2,b12,b1) (b12+b2)*(2*x-1)+b1
W_a = function(x,b2,b12,b1) (b12-b2)*(2*x-1)-b1
W_m = function(x,b2,b12,b1) 2*b2*x*(2*x-1)+(b12-b2)*(2*x-1)-b1+2*b1*x

wplot = function(b2,b12,b1,pch) {
  p = ggplot(NULL,aes(x=2,y=2)) + 
    geom_point(pch=pch,size=3) + ylim(0,2) + xlim(0,1) + theme_classic() +
    geom_function(fun=function(x,b2,b12,b1) (b12+b2)*(2*x-1)+b1+1, col=grey(0.0,1.0),args=list(b2,b12,b1)) +
    geom_function(fun=function(x,b2,b12,b1) (b12-b2)*(2*x-1)-b1+1, col=grey(0.0,0.5),args=list(b2,b12,b1)) + 
    geom_function(fun=function(x,b2,b12,b1) 2*b2*x*(2*x-1)+(b12-b2)*(2*x-1)-b1+2*b1*x+1, lty=2,args=list(b2,b12,b1)) + 
    ylab("Absolute fitness") + xlab("Frequency of AA genotypes") +
    labs(title=bquote(beta[s]*"="*.(b1)*"; "*beta[n]*"="*.(b12)*"; "*beta[sxn]*"="*.(b2)))
  return(p)
} 

pa = wplot(b2=-0.3,b12=0,b1=0,pch=16)
pa = pa + geom_text(aes(x=0.025,y=2.0),label="(black): AA genotypes",size=3,hjust=0) + 
  geom_text(aes(x=0.025,y=1.8),label="(gray): aa genotypes",colour=grey(0.5,1),size=3,hjust=0) +
  geom_text(aes(x=0.025,y=1.6),label="(dashed): mean",size=3,hjust=0)
pb = wplot(b2=-0.3,b12=0,b1=0.2,pch=16)
pc = wplot(b2=-0.3,b12=0.2,b1=0,pch=16)
pd = wplot(b2=-0.3,b12=0.2,b1=0.2,pch=16)
pe = wplot(b2=0.3,b12=0,b1=0,pch=16)
pf = wplot(b2=0.3,b12=0,b1=0.2,pch=16)
pg = wplot(b2=0.3,b12=0.2,b1=0,pch=16)
ph = wplot(b2=0.3,b12=0.2,b1=0.2,pch=16)

p = (pa | pb | pc | pd) / (pe | pf | pg | ph) + plot_annotation(tag_levels = "a")
ggsave(p,filename="../figs/FDS.pdf",width=11,height=5.5)

###########
# Figure 2 & Figure S1
ggMan = function(chr,pos,p) {
  chr_rep = table(chr)
  alphas = 1-(chr/2 - chr %/% 2)
  cols = rgb(0,0,alphas,alphas)
  x = coord(chr,pos)
  y = -log10(p)
  man = ggplot(NULL,aes(x=x$coord,y=y)) + geom_point(colour=cols) + theme_classic() + 
    scale_x_continuous(name="Chromosomes", breaks=x$tic, labels=names(chr_rep)) +
    ylab(expression(-log[10]*(italic(p)))) + geom_hline(yintercept=-log10(0.05/length(p)),lty=2,colour="black")
  return(man)
}

ggQQ = function(chr,pos,p) {
  chr_rep = table(chr)
  alphas = 1-(chr/2 - chr %/% 2)
  cols = rgb(0,0,alphas,alphas)
  cols = cols[order(p,decreasing=FALSE)]
  
  o = -log(sort(p,decreasing=FALSE),10)
  e = -log(ppoints(length(p)),10)
  
  qq = ggplot(data=NULL, mapping=aes(x=e,y=o))+
    geom_point(colour=cols)+
    geom_abline(intercept=0,slope=1,linetype="dashed")+
    theme_classic()+
    xlab(expression("Expected "*-log[10](p)))+ylab(expression("Observed "*-log[10](p)))
  return(qq)
}

ggEHH = function(chr,pos,p) {
  chr_rep = table(chr)
  alphas = 1-(chr/2 - chr %/% 2)
  cols = rgb(0,0,alphas,alphas)
  x = coord(chr,pos)
  y = p
  man = ggplot(NULL,aes(x=x$coord,y=y)) + geom_point(colour=cols) + theme_classic() + 
    scale_x_continuous(name="Chromosomes", breaks=x$tic, labels=names(chr_rep)) +
    geom_hline(yintercept=quantile(y,0.975,na.rm=TRUE),lty=2,colour="black")
  return(man)
}

perm_out = readRDS(file="../output/Wuest2022perm/NeiGWAS_Wuest_et_al_sim_all_biomassZ_perm.rds")
hist(-log10(perm_out[,3]))
abline(v=quantile(-log10(perm_out[,3]),0.95),lty=2)
quantile(-log10(perm_out[,3]),0.95)
permp = ggplot(NULL,aes(x=-log10(perm_out[,3]))) + geom_histogram() + 
  geom_vline(xintercept=quantile(-log10(perm_out[,3]),0.95),lty=2) + 
  xlab(expression(-log[10]*(italic(p)))) + theme_classic()

d1 = readRDS(file="../output/NeiGWAS_Wuest_et_al_sim_all_biomassZ.rds")
man11 = ggMan(d1$chr,d1$pos,d1$p_s)
qq11 = ggQQ(d1$chr,d1$pos,d1$p_s)

man12 = ggMan(d1$chr,d1$pos,d1$p_n)
qq12 = ggQQ(d1$chr,d1$pos,d1$p_n)

man13 = ggMan(d1$chr,d1$pos,d1$p_sim) 
man13 = man13 + geom_hline(yintercept=quantile(-log10(perm_out[,3]),0.95),lty=2,colour="grey")
qq13 = ggQQ(d1$chr,d1$pos,d1$p_sim)

d2 = readRDS(file="../output/NeiGWAS_Wuest_et_al_sim_all_rosette.rds")
man21 = ggMan(d2$chr,d2$pos,d2$p_s) 
qq21 = ggQQ(d2$chr,d2$pos,d2$p_s)

man22 = ggMan(d2$chr,d2$pos,d2$p_n)
qq22 = ggQQ(d2$chr,d2$pos,d2$p_n)

man23 = ggMan(d2$chr,d2$pos,d2$p_sim)
qq23 = ggQQ(d2$chr,d2$pos,d2$p_sim)

d3 = readRDS(file="../output/NeiGWAS_Wuest_et_al_sim_all_flowering.rds")
man31 = ggMan(d3$chr,d3$pos,d3$p_s) 
qq31 = ggQQ(d3$chr,d3$pos,d3$p_s)

man32 = ggMan(d3$chr,d3$pos,d3$p_n) 
qq32 = ggQQ(d3$chr,d3$pos,d3$p_n)

man33 = ggMan(d3$chr,d3$pos,d3$p_sim) 
qq33 = ggQQ(d3$chr,d3$pos,d3$p_sim)

man = (man11 + man21 + man31) / (man12 + man22 + man23) / (man13 + man23 + man33)
ggsave(man,filename="../figs/wuest_manhattan.jpg",dpi=600,width=8,height=4)

qq = ((qq11+labs(title="biomass",subtitle="self")) + (qq21+labs(title="flowering time",subtitle="self")) + (qq31+labs(title="rosette size",subtitle="self"))) / ((qq12+labs(subtitle="neighbor")) + (qq22+labs(subtitle="neighbor")) + (qq23+labs(subtitle="neighbor"))) / ((qq13+labs(subtitle="self x nei.")) + (qq23+labs(subtitle="self x nei.")) + (qq33)+labs(subtitle="self x nei.")) + plot_annotation(tag_levels = "a")
ggsave(qq,filename="../figs/wuest_qq.jpg",dpi=600,width=9,height=9)

# perm out ver2
pvec = c()
for(i in 1:199) {
  out = readRDS(paste0("../output/Wuest2022perm/grotate/NeiGWAS_Wuest_et_al_sim_all_biomassZ_",i,".rds"))
  pvec = c(pvec,min(out$p_sim))
}
pvec = -log10(pvec)
hist(pvec)
abline(v=quantile(pvec,0.95),lty=2)
rotp = ggplot(NULL,aes(x=pvec)) + geom_histogram() + 
  geom_vline(xintercept=quantile(pvec,0.95),lty=2) + 
  xlab(expression(-log[10]*(italic(p)))) + theme_classic()

permh = permp + rotp + plot_annotation(tag_levels="a")
ggsave(permh,filename="Wuest_perm.pdf",width=6,height=3)

################
# Figure 2b & Figure SX
th = quantile(perm_out[,3],0.05)
d = readRDS(file="../output/NeiGWAS_Wuest_et_al_sim_all_biomassZ_beta.rds")
d[d$p_sim<th,]
write.csv(d[d$p_sim<th,],"../output/betasxn_sig.csv")

d[order(d$p_sim),]
d[order(d$p_sim),][1,]
plot(jitter(g_sim[,165819]),pheno2$biomass)
plot(factor(g_sim[,165819]),pheno2$biomass)
plot(pheno2$biomass~factor(factor(g_self[,165819]):factor(g_nei[,165819])))
aggregate(pheno2$biomass~factor(factor(g_self[,165819]):factor(g_nei[,165819])),data=NULL,mean)


resid = lm(pheno2$biomass~g_self[,165819])$residuals
t.test(pheno2$biomass~factor(g_sim[,165819]))
t.test(resid~factor(g_sim[,165819]))
plot(resid~factor(g_sim[,165819]))

# chr5:2838468 (Montazeaud et al. 2023)
plot(jitter(g_sim[,165823]),pheno2$biomass)
plot(factor(g_sim[,165823]),pheno2$biomass)
plot(pheno2$biomass~factor(factor(g_self[,165823]):factor(g_nei[,165823])))


resid = lm(pheno2$biomass~g_self[,165835]+g_nei[,165835])$residuals
t.test(resid~factor(g_sim[,165835]))
t.test(pheno2$biomass~factor(g_self[,11063]))
t.test(pheno2$biomass~factor(g_self[,1473]))

m2 = mean(pheno2$biomass[g_sim[,165819]==1])
sd2 = sd(pheno2$biomass[g_sim[,165819]==1]) #/sqrt(length(pheno2$biomass[g_sim[,165819]==1]))

m1 = mean(pheno2$biomass[g_sim[,165819]==-1])
sd1 = sd(pheno2$biomass[g_sim[,165819]==-1]) #/sqrt(length(pheno2$biomass[g_sim[,165819]==-1]))

data = data.frame(pheno2$biomass,factor(g_sim[,165819]))
colnames(data) = c("biomass","similarity")
p1 = ggplot(data=data,mapping=aes(x=similarity,y=biomass)) + geom_jitter(col=grey(0.5,0.25)) + 
  geom_violin(alpha=0.5) + geom_boxplot(width=0.3,outlier.shape=NA,alpha=0.5) + theme_classic()
  #geom_point(aes(x=c(-1,1),y=c(m1,m2))) + 
  #geom_errorbar(aes(x=factor(c(-1,1)),ymin=c(m1-sd1,m2-sd2),y=c(m1+sd1,m2+sd2))) + theme_classic()

plot(pheno2$biomass~factor(factor(g_self[,165819]):factor(g_nei[,165819])))
aggregate(pheno2$biomass~factor(factor(g_self[,165819]):factor(g_nei[,165819])),data=NULL,mean)

data = data.frame(pheno2$biomass,factor(factor(g_self[,165819]):factor(g_nei[,165819])))
colnames(data) = c("biomass","self_neig")
p2 = ggplot(data=data,mapping=aes(x=self_neig,y=biomass)) + geom_jitter(col=grey(0.5,0.25)) + 
  geom_violin(alpha=0.5) + geom_boxplot(width=0.3,outlier.shape=NA,alpha=0.5) + theme_classic()

f = (g_self[,165819]+g_nei[,165819])/4+0.5
lm(pheno2[which(g_self[,165819]==1),"biomass"]~f[which(g_self[,165819]==1)])
lm(pheno2[which(g_self[,165819]==-1),"biomass"]~f[which(g_self[,165819]==-1)])

plot(jitter(f),pheno2$biomass,pch=16,col=grey(0.5,0.5),las=1)
abline(a=539.17,b=-40.85)
abline(a=505.31,b=39.27,lty=2)

tar = 165823 # topSNP: 165819; Germain's SNP: 165823; Chr1 SNP: 1473, 11063 looks imcomplete but fine, 11074 11071 11086 seem selection effects
f = (g_self[,tar]+g_nei[,tar])/4+0.5
plot(jitter(f),pheno2$biomass,pch=g_self[,tar]+2,col=grey(0.5,0.5),las=1)
res1 = lm(pheno2[which(g_self[,tar]==1),"biomass"]~f[which(g_self[,tar]==1)])
res2 = lm(pheno2[which(g_self[,tar]==-1),"biomass"]~f[which(g_self[,tar]==-1)])
abline(a=res1$coefficients[1],b=res1$coefficients[2])
abline(a=res2$coefficients[1],b=res2$coefficients[2],lty=2)

p1 = ggplot(data=NULL,mapping=aes(x=jitter(f),y=pheno2$biomass)) + geom_point(pch=g_self[,tar]+2,col="grey") + 
  geom_abline(intercept=res1$coefficients[1],slope=res1$coefficients[2]) +
  geom_abline(intercept=res2$coefficients[1],slope=res2$coefficients[2],lty=2) + 
  ylab("biomass (mg)") + xlab("freq. of alternative alleles") + 
  geom_text(aes(x=0.7,y=1100),label="(dashed): Ref.",size=2,hjust=0) + 
  geom_text(aes(x=0.7,y=1050),label="(solid): Alt.",size=2,hjust=0) + 
  theme_classic()
  
data = data.frame(pheno2$biomass,factor(factor(g_self[,tar]):factor(g_nei[,tar])))
colnames(data) = c("biomass","self_neig")
p2 = ggplot(data=data,mapping=aes(x=self_neig,y=biomass)) + geom_jitter(pch=g_self[,tar]+2,col=grey(0.5,0.25)) + 
  geom_violin(alpha=0.5) + geom_boxplot(width=0.3,outlier.shape=NA,alpha=0.5) + 
  scale_x_discrete(labels=c("ref:ref","ref:alt","alt:ref","alt:alt")) +
  ylab("biomass (mg)") + xlab("self:neighbor") + 
  geom_point(aes(x=3.5,y=1100),pch=1) + 
  geom_point(aes(x=3.5,y=1050),pch=3) + 
  geom_text(aes(x=3.65,y=1100),label="Ref.",size=2.5,hjust=0) + 
  geom_text(aes(x=3.65,y=1050),label="Alt.",size=2.5,hjust=0) + 
  theme_classic() 
  

p3 = ggplot(data=NULL,mapping=aes(x=factor(g_sim[,tar]),y=pheno2$biomass)) + geom_jitter(pch=g_self[,tar]+2,col=grey(0.5,0.25)) + 
  geom_violin(alpha=0.5) + geom_boxplot(width=0.3,outlier.shape=NA,alpha=0.5) + ylab("biomass (mg)") + xlab("") + 
  scale_x_discrete(labels=c("mixture","monoculture")) + theme_classic()
  # theme(legend.background = element_rect(fill = NA, colour = NA),
  #       plot.background = element_rect(fill=NA, color=NA),
  #       panel.background = element_rect(fill = NA, colour = "black"))  


plot(jitter(g_sim[,tar]),pheno2$biomass)
plot(factor(g_sim[,tar]),pheno2$biomass)
plot(pheno2$biomass~factor(factor(g_self[,tar]):factor(g_nei[,tar])))

ggsave(p2+p1+p3,filename="../figs/wuest_boxplot2.pdf",width=7,height=3)

#########################
# Figure 3a: output of fine mapping
d2 = readRDS("../output/NeiGWAS_Wuest_et_al_sim_all_biomassZ_impMAF5.rds")
x = data.frame(d2$chr,d2$pos,d2$p_sim)
colnames(x) = c("chr","pos","p")

d2[order(d2$p_sim),]
d2[d2$chr==5&d2$pos==2832310,]
plot(d2[c((1499421-500):(1499421+1200)),"pos"],-log10(d2[c((1499421-500):(1499421+1200)),"p_sim"]),
     ylab="-log10(p)",xlab="position (bp)")

mnchr5 = ggplot(NULL, aes(x=d2[c((1499421-500):(1499421+1000)),"pos"],y=-log10(d2[c((1499421-500):(1499421+1000)),"p_sim"]))) + 
  geom_point(colour="blue") + ylab(expression(-log[10]*(italic(p)))) + xlab("Position (bp)") + theme_bw()

chr5p = (mnchr5 / plot_spacer() / betap / ehhp)
ggsave(chr5p, filename="Wuest_chr5.pdf",width=7,height=6)

##############
# Figure 3c,d & Figure SX: scan output
# load and merge scan results
chr1_hh = read.csv("../output/scan_ihs98_chr1.csv")
chr2_hh = read.csv("../output/scan_ihs98_chr2.csv")
chr3_hh = read.csv("../output/scan_ihs98_chr3.csv")
chr4_hh = read.csv("../output/scan_ihs98_chr4.csv")
chr5_hh = read.csv("../output/scan_ihs98_chr5.csv")

ehh_all = rbind(chr1_hh,chr2_hh,chr3_hh,chr4_hh,chr5_hh)

ehh_all = na.omit(ehh_all)
ehh_out = ehh_all[ehh_all$IHS>quantile(ehh_all$IHS,0.95),]


chr1_beta = read.table("../output/BetaScanChr1out.txt", header=TRUE)
chr2_beta = read.table("../output/BetaScanChr2out.txt", header=TRUE)
chr3_beta = read.table("../output/BetaScanChr3out.txt", header=TRUE)
chr4_beta = read.table("../output/BetaScanChr4out.txt", header=TRUE)
chr5_beta = read.table("../output/BetaScanChr5out.txt", header=TRUE)
Chr = c(rep(1,nrow(chr1_beta)), rep(2,nrow(chr2_beta)), rep(3,nrow(chr3_beta)), rep(4,nrow(chr4_beta)), rep(5,nrow(chr5_beta)))

beta_all = rbind(chr1_beta,chr2_beta,chr3_beta,chr4_beta,chr5_beta)
beta_all = data.frame(Chr,beta_all)
beta_out = beta_all[beta_all$Beta1>quantile(beta_all$Beta1,0.95,na.rm=TRUE),]


top = d[-log10(d$p_sim)>7.993031,]
paste0(top$chr,"-",top$pos)
intersect(paste0(beta_out$Chr,"-",beta_out$Position),paste0(ehh_out$CHR,"-",ehh_out$POSITION))[2000:2626]


beta_out[beta_out$Chr==5,]
ehh_out[ehh_out$CHR==5,][500:4500,]

beta_all[(beta_all$Chr==5&beta_all$Position>2800000)&beta_all$Chr==5&beta_all$Position<2940000,]

plot(beta_all[(beta_all$Chr==5&beta_all$Position>2800000)&beta_all$Chr==5&beta_all$Position<2900000,"Position"],
     beta_all[(beta_all$Chr==5&beta_all$Position>2800000)&beta_all$Chr==5&beta_all$Position<2900000,"Beta1"],
     type="l",ylab="BETA",xlab="position (bp)")
abline(h=quantile(beta_all$Beta1,0.975,na.rm=TRUE),lty=2)

betap = ggplot(NULL, aes(x=beta_all[(beta_all$Chr==5&beta_all$Position>2800000)&beta_all$Chr==5&beta_all$Position<2900000,"Position"],
                         y=beta_all[(beta_all$Chr==5&beta_all$Position>2800000)&beta_all$Chr==5&beta_all$Position<2900000,"Beta1"])) + 
  geom_line() + ylab("BETA(1)") + xlab("Position (bp)") + theme_bw() +
  geom_hline(yintercept=quantile(beta_all$Beta1,0.975,na.rm=TRUE),lty=2)


plot(ehh_all[(ehh_all$CHR==5&ehh_all$POSITION>2800000)&ehh_all$CHR==5&ehh_all$POSITION<2900000,"POSITION"],
     ehh_all[(ehh_all$CHR==5&ehh_all$POSITION>2800000)&ehh_all$CHR==5&ehh_all$POSITION<2900000,"IHS"],
     type="l",ylab="iHS",xlab="position (bp)")
abline(h=quantile(ehh_all$IHS,0.975,na.rm=TRUE),lty=2)

ehhp = ggplot(NULL, aes(x=ehh_all[(ehh_all$CHR==5&ehh_all$POSITION>2800000)&ehh_all$CHR==5&ehh_all$POSITION<2900000,"POSITION"],
                        y=ehh_all[(ehh_all$CHR==5&ehh_all$POSITION>2800000)&ehh_all$CHR==5&ehh_all$POSITION<2900000,"IHS"])) + 
  geom_line() + ylab("iHS") + xlab("Position (bp)") + theme_bw() +
  geom_hline(yintercept=quantile(ehh_all$IHS,0.975,na.rm=TRUE),lty=2)


hbeta = ggplot(beta_all,aes(x=Beta1)) + geom_histogram() + ylab("BETA(1)") +
  geom_vline(xintercept=quantile(beta_all$Beta1,0.975,na.rm=TRUE),lty=2) + 
  coord_flip() + theme_classic()

hehh = ggplot(ehh_all,aes(x=IHS)) + geom_histogram() + ylab("iHS") +
  geom_vline(xintercept=quantile(ehh_all$IHS,0.975,na.rm=TRUE),lty=2) + 
  coord_flip() + theme_classic()

mehh = ggEHH(chr=ehh_all$CHR,pos=ehh_all$POSITION,p=ehh_all$IHS) + ylab("iHS")
mbeta = ggEHH(chr=beta_all$Chr,pos=beta_all$Position,p=beta_all$Beta1) + ylab("BETA(1)")

scanp = (hbeta + mbeta + plot_layout(widths=c(1,3))) / (hehh + mehh + plot_layout(widths=c(1,3)))
scanp = scanp + plot_annotation(tag_levels = "a")
ggsave(scanp,filename="Wuest_scan.jpg",width=8,height=4,dpi=300)


##############
# Figure SX: haplotype check

# ten tester genotypes
testers = names(table(pheno$AccID_Tester_posA)[table(pheno$AccID_Tester_posA)>20])
# 6899 6906 6909 6911 6932 6962 7328 7378 8387 8412 
# 190  193  192  189  193  193  191  191  195  196 

par(mfcol=c(3,4))
for(i in testers) {
  tar = which(pheno2$focal==paste0("X",i))
  pval = round(t.test(pheno2[tar,]$biomass~factor(g_sim[tar,165819]))$p.value,2)
  plot(pheno2[tar,]$biomass~factor(g_sim[tar,165819]),main=paste(i,pval))
}

image(g_self[paste0("X",testers),c((165819-100):(165819+100))])

g_self[paste0("X",testers),165819]
# 6909: Col-0

# imputed genotypes
gimp = readRDS("../data/SEW2022_sub_snpsMAF5.rds")
posimp = readRDS("../data/SEW_sub_posMAF5.rds")

posimp[(posimp$chr==5&posimp$pos>2800000)&(posimp$chr==5&posimp$pos<2900000),]
whichg = which((posimp$chr==5&posimp$pos>2800000)&(posimp$chr==5&posimp$pos<2900000))

png("../figs/WuestHaplo.png",width=6,height=6,res=220,units="in")
par(bg=NA)
image(as.matrix(gimp[whichg,]),ylab="accessions",xlab="position (bp)",xaxt="n",yaxt="n")
dev.off()

ld = cor(t(as.matrix(gimp[whichg,])))
png("../figs/WuestLD.png",width=6,height=6,res=220,units="in")
par(bg=NA)
heatmap(x=ld,xlab="position (bp)",ylab="position (bp)")
dev.off()

posimp[(posimp$chr==5&posimp$pos>2800000)&(posimp$chr==5&posimp$pos<2900000),]

res = prcomp(ld)
png("../figs/WuestHaploLDpca.png",width=12,height=12,res=220,units="in")
par(mfcol=c(2,2),bg=NA)
image(as.matrix(gimp[whichg,]),ylab="accessions",xlab="position (bp)",xaxt="n",yaxt="n",main="(a) haplotypes")
image(x=posimp[whichg,"pos"],y=posimp[whichg,"pos"],ld,xlab="position (bp)",ylab="position (bp)",main="(b) LD")
plot(NULL)
plot(summary(res)$importance[3,1:6],ylim=c(0,1),las=1,type="b", ylab="Cumulative prop. of variance explained",xlab="Principal component",main="(c) PCA")
dev.off()

# png("../figs/WuestPCA.png",width=5,height=5,res=220,units="in")
# par(bg=NA)
# plot(summary(res)$importance[3,1:6],ylim=c(0,1),las=1,type="b",
#      ylab="Cumulative prop. of variance explained",
#      xlab="PCA")
# # points(summary(res)$importance[2,1:6],pch=20)
# dev.off()

###############
# Genome inflation factor
GIF = function(d,type="p_sim") {
  numer = median(-log10(d[,type]),0.5)
  denom = median(-log10(ppoints(length(d[,type]))))
  lambda_gc = numer/denom
  return(lambda_gc)
}

mgif= matrix(c(
GIF(d=readRDS(file="../output/NeiGWAS_Wuest_et_al_sim_all_biomassZ.rds"),type="p_s"),
GIF(d=readRDS(file="../output/NeiGWAS_Wuest_et_al_sim_all_flowering.rds"),type="p_s"),
GIF(d=readRDS(file="../output/NeiGWAS_Wuest_et_al_sim_all_rosette.rds"),type="p_s"),

GIF(d=readRDS(file="../output/NeiGWAS_Wuest_et_al_sim_all_biomassZ.rds"),type="p_n"),
GIF(d=readRDS(file="../output/NeiGWAS_Wuest_et_al_sim_all_flowering.rds"),type="p_n"),
GIF(d=readRDS(file="../output/NeiGWAS_Wuest_et_al_sim_all_rosette.rds"),type="p_n"),

GIF(d=readRDS(file="../output/NeiGWAS_Wuest_et_al_sim_all_biomassZ.rds")),
GIF(d=readRDS(file="../output/NeiGWAS_Wuest_et_al_sim_all_flowering.rds")),
GIF(d=readRDS(file="../output/NeiGWAS_Wuest_et_al_sim_all_rosette.rds"))
), 3,3)

mgif = as.data.frame(mgif)

colnames(mgif) = c("self","nei","selfxnei")
rownames(mgif) = c("biomass","flowering","rosette")

write.csv(mgif,"../output/Wuest_GIF.csv")

