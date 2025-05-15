library(tidyverse)
library(patchwork)
library(rNeighborGWAS)
library(gaston)
source("./script/coord.R")

####################
# reload genotypes and phenotypes
pheno = read.csv("./pheno/competition.csv")

pheno = filter(pheno,CommunityType!="single")
pheno2 = gather(pheno,"Biomass_mg_posA","Biomass_mg_posB",key="posAB",value="biomass")

focal = c(pheno$AccID_Tester_posA,pheno$AccID_Target_posB)
neighbor = c(pheno$AccID_Target_posB,pheno$AccID_Tester_posA)
pheno2 = data.frame(focal,neighbor,pheno2)

pheno2$focal = paste0("X",pheno2$focal)
pheno2$neighbor = paste0("X",pheno2$neighbor)
pheno2 = filter(pheno2, pheno2$focal!="XNA"|pheno2$neighbor!="XNA")
pheno2 = pheno2[-which(is.na(pheno2$biomass)),]

#load genotype data
g = read.csv("./geno/call_method_75_TAIR9_250k.csv.gz",header=TRUE,skip=1)
line_names = colnames(g)[-c(1:2)]
info = g[,1:2]
g = g[,-c(1:2)]

#convert rare alleles into -1 degit
snp_degit = function(vec) {
  if(table(vec)[1]>table(vec)[2]) {
    vec[which(vec==levels(factor(vec))[1])] = 0
    vec[which(vec==levels(factor(vec))[2])] = 1
  } else {
    vec[which(vec==levels(factor(vec))[1])] = 1
    vec[which(vec==levels(factor(vec))[2])] = 0
  }
  return(as.numeric(vec))
}

g_bin = apply(g,1,snp_degit)
g_bin = t(g_bin)
colnames(g_bin) = line_names

geno = list(pos=info,g=g_bin)

rm(g_bin); rm(g)
gc();gc()

geno$g = geno$g[,levels(factor(c(pheno2$focal,pheno2$neighbor)))]

# MAF cutoff 5%
AF = apply(geno$g,1,sum)/ncol(geno$g)
TF = which(AF>0.05&AF<0.95)
AF = AF[which(AF>0.05&AF<0.95)]

geno$g = geno$g[TF,]
geno$g[geno$g==0] = -1

geno$pos = geno$pos[TF,]

g_self = geno$g[,pheno2$focal]
g_nei = geno$g[,pheno2$neighbor]
g_sim = mapply(function(x) { return(geno$g[,pheno2$focal[x]]*geno$g[,pheno2$neighbor[x]]) }, 1:nrow(pheno2))

g_self = t(g_self)
g_nei = t(g_nei)
g_sim = t(g_sim)


############
# main figure: numerical example
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
ggsave(p,filename="FDS.pdf",width=11,height=5.5)

###########
# Manhattan & QQ plot
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

# p-values of permutation tests
# p_neig: neighbor randomization; p_grot: genome rotation
perm_out = read.csv("./output/p_perm.csv")
permp = ggplot(NULL,aes(x=-log10(perm_out$p_neig))) + geom_histogram() + 
  geom_vline(xintercept=quantile(-log10(perm_out$p_neig),0.95),lty=2) + 
  xlab(expression(-log[10]*(italic(p)))) + theme_classic()

d1 = readRDS(file="./output/NeiGWAS_Wuest_et_al_sim_all_biomassZ.rds")
man11 = ggMan(d1$chr,d1$pos,d1$p_s)
qq11 = ggQQ(d1$chr,d1$pos,d1$p_s)

man12 = ggMan(d1$chr,d1$pos,d1$p_n)
qq12 = ggQQ(d1$chr,d1$pos,d1$p_n)

man13 = ggMan(d1$chr,d1$pos,d1$p_sim) 
man13 = man13 + geom_hline(yintercept=quantile(-log10(perm_out$p_neig),0.95),lty=2,colour="grey")
qq13 = ggQQ(d1$chr,d1$pos,d1$p_sim)

d2 = readRDS(file="./output/NeiGWAS_Wuest_et_al_sim_all_rosetteZ.rds")
man21 = ggMan(d2$chr,d2$pos,d2$p_s) 
qq21 = ggQQ(d2$chr,d2$pos,d2$p_s)

man22 = ggMan(d2$chr,d2$pos,d2$p_n)
qq22 = ggQQ(d2$chr,d2$pos,d2$p_n)

man23 = ggMan(d2$chr,d2$pos,d2$p_sim)
qq23 = ggQQ(d2$chr,d2$pos,d2$p_sim)

d3 = readRDS(file="./output/NeiGWAS_Wuest_et_al_sim_all_floweringZ.rds")
man31 = ggMan(d3$chr,d3$pos,d3$p_s) 
qq31 = ggQQ(d3$chr,d3$pos,d3$p_s)

man32 = ggMan(d3$chr,d3$pos,d3$p_n) 
qq32 = ggQQ(d3$chr,d3$pos,d3$p_n)

man33 = ggMan(d3$chr,d3$pos,d3$p_sim) 
qq33 = ggQQ(d3$chr,d3$pos,d3$p_sim)

man = (man11 + man21 + man31) / (man12 + man22 + man23) / (man13 + man23 + man33)
ggsave(man,filename="wuest_manhattan.jpg",dpi=600,width=8,height=4)

qq = ((qq11+labs(title="biomass",subtitle="self")) + (qq21+labs(title="flowering time",subtitle="self")) + (qq31+labs(title="rosette size",subtitle="self"))) / ((qq12+labs(subtitle="neighbor")) + (qq22+labs(subtitle="neighbor")) + (qq23+labs(subtitle="neighbor"))) / ((qq13+labs(subtitle="self x nei.")) + (qq23+labs(subtitle="self x nei.")) + (qq33)+labs(subtitle="self x nei.")) + plot_annotation(tag_levels = "a")
ggsave(qq,filename="wuest_qq.jpg",dpi=600,width=9,height=9)

# permutation histogram for the genome rotation scheme
rotp = ggplot(NULL,aes(x=-log10(perm_out$p_grot))) + geom_histogram() + 
  geom_vline(xintercept=quantile(-log10(perm_out$p_grot),0.95),lty=2) + 
  xlab(expression(-log[10]*(italic(p)))) + theme_classic()

permh = permp + rotp + plot_annotation(tag_levels="a")
ggsave(permh,filename="Wuest_perm.pdf",width=6,height=3)

################
# Boxplots for main and supplementary figures

# check significant SNPs
d = readRDS(file="./output/NeiGWAS_Wuest_et_al_sim_all_biomassZ.rds")
th = quantile(perm_out$p_neig,0.05) # set threshold
d[d$p_sim<th,]

# top SNP on chr5: '165819' for a main figure
# top SNP on chr1: '11063' for supplementary figures
# Montazeaud et al.'s (2023) SNP: '165823' for supplementary figures
tar = 165819 # select a target SNP
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
  scale_x_discrete(labels=c("biallelic","monoallelic")) + theme_classic()

ggsave(p2+p1+p3,filename="wuest_boxplot1.pdf",width=7,height=3)

#########################
# Main figure: fine mapping
d2 = readRDS("./output/NeiGWAS_Wuest_et_al_sim_all_biomassZ_impMAF5.rds")
x = data.frame(d2$chr,d2$pos,d2$p_sim)
colnames(x) = c("chr","pos","p")

# check the position near the most significant SNP
d2[order(d2$p_sim),]
d2[d2$chr==5&d2$pos==2832310,]

mnchr5 = ggplot(NULL, aes(x=d2[c((1499421-500):(1499421+1000)),"pos"],y=-log10(d2[c((1499421-500):(1499421+1000)),"p_sim"]))) + 
  geom_point(colour="blue") + ylab(expression(-log[10]*(italic(p)))) + xlab("Position (bp)") + theme_bw()

chr5p = (mnchr5 / plot_spacer() / betap / ehhp)
ggsave(chr5p, filename="Wuest_chr5.pdf",width=7,height=6)

##############
# genome scan for main and supplementary figures
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

# load and merge scan results
chr1_hh = read.csv("./output/scan_ihs98_chr1.csv.gz")
chr2_hh = read.csv("./output/scan_ihs98_chr2.csv.gz")
chr3_hh = read.csv("./output/scan_ihs98_chr3.csv.gz")
chr4_hh = read.csv("./output/scan_ihs98_chr4.csv.gz")
chr5_hh = read.csv("./output/scan_ihs98_chr5.csv.gz")

ehh_all = rbind(chr1_hh,chr2_hh,chr3_hh,chr4_hh,chr5_hh)

ehh_all = na.omit(ehh_all)
ehh_out = ehh_all[ehh_all$IHS>quantile(ehh_all$IHS,0.95),]

chr1_beta = read.table("./output/BetaScanChr1out.txt", header=TRUE)
chr2_beta = read.table("./output/BetaScanChr2out.txt", header=TRUE)
chr3_beta = read.table("./output/BetaScanChr3out.txt", header=TRUE)
chr4_beta = read.table("./output/BetaScanChr4out.txt", header=TRUE)
chr5_beta = read.table("./output/BetaScanChr5out.txt", header=TRUE)
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

betap = ggplot(NULL, aes(x=beta_all[(beta_all$Chr==5&beta_all$Position>2800000)&beta_all$Chr==5&beta_all$Position<2900000,"Position"],
                         y=beta_all[(beta_all$Chr==5&beta_all$Position>2800000)&beta_all$Chr==5&beta_all$Position<2900000,"Beta1"])) + 
  geom_line() + ylab("BETA(1)") + xlab("Position (bp)") + theme_bw() +
  geom_hline(yintercept=quantile(beta_all$Beta1,0.975,na.rm=TRUE),lty=2)


ehhp = ggplot(NULL, aes(x=ehh_all[(ehh_all$CHR==5&ehh_all$POSITION>2800000)&ehh_all$CHR==5&ehh_all$POSITION<2900000,"POSITION"],
                        y=ehh_all[(ehh_all$CHR==5&ehh_all$POSITION>2800000)&ehh_all$CHR==5&ehh_all$POSITION<2900000,"IHS"])) + 
  geom_line() + ylab("iHS") + xlab("Position (bp)") + theme_bw() +
  geom_hline(yintercept=quantile(ehh_all$IHS,0.975,na.rm=TRUE),lty=2)

hbeta = ggplot(beta_all,aes(x=Beta1)) + geom_histogram() + xlab("BETA(1)") +
  geom_vline(xintercept=quantile(beta_all$Beta1,0.975,na.rm=TRUE),lty=2) + 
  coord_flip() + theme_classic()

hehh = ggplot(ehh_all,aes(x=IHS)) + geom_histogram() + xlab("iHS") +
  geom_vline(xintercept=quantile(ehh_all$IHS,0.975,na.rm=TRUE),lty=2) + 
  coord_flip() + theme_classic()

mehh = ggEHH(chr=ehh_all$CHR,pos=ehh_all$POSITION,p=ehh_all$IHS) + ylab("iHS")
mbeta = ggEHH(chr=beta_all$Chr,pos=beta_all$Position,p=beta_all$Beta1) + ylab("BETA(1)")

scanp = (hbeta + mbeta + plot_layout(widths=c(1,3))) / (hehh + mehh + plot_layout(widths=c(1,3)))
scanp = scanp + plot_annotation(tag_levels = "a")
ggsave(scanp,filename="Wuest_scan.jpg",width=8,height=4,dpi=300)


##############
# haplotype, LD, and PCA for supplementary figures

# load imputed genotypes
gimp = readRDS("./geno/SEW2022_sub_snpsMAF5.rds")
posimp = readRDS("./geno/SEW_sub_posMAF5.rds")

posimp[(posimp$chr==5&posimp$pos>2800000)&(posimp$chr==5&posimp$pos<2900000),]
whichg = which((posimp$chr==5&posimp$pos>2800000)&(posimp$chr==5&posimp$pos<2900000))

posimp[(posimp$chr==5&posimp$pos>2800000)&(posimp$chr==5&posimp$pos<2900000),]

ld = cor(t(as.matrix(gimp[whichg,])))
res = prcomp(ld)
png("WuestHaploLDpca.png",width=12,height=12,res=220,units="in")
par(mfcol=c(2,2),bg=NA)
image(as.matrix(gimp[whichg,]),ylab="accessions",xlab="position (bp)",xaxt="n",yaxt="n",main="(a) haplotypes")
image(x=posimp[whichg,"pos"],y=posimp[whichg,"pos"],ld,xlab="position (bp)",ylab="position (bp)",main="(b) LD")
try(plot(NULL))
plot(summary(res)$importance[3,1:6],ylim=c(0,1),las=1,type="b", ylab="Cumulative prop. of variance explained",xlab="Principal component",main="(c) PCA")
dev.off()


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
GIF(d=readRDS(file="../output/NeiGWAS_Wuest_et_al_sim_all_floweringZ.rds"),type="p_s"),
GIF(d=readRDS(file="../output/NeiGWAS_Wuest_et_al_sim_all_rosetteZ.rds"),type="p_s"),

GIF(d=readRDS(file="../output/NeiGWAS_Wuest_et_al_sim_all_biomassZ.rds"),type="p_n"),
GIF(d=readRDS(file="../output/NeiGWAS_Wuest_et_al_sim_all_floweringZ.rds"),type="p_n"),
GIF(d=readRDS(file="../output/NeiGWAS_Wuest_et_al_sim_all_rosetteZ.rds"),type="p_n"),

GIF(d=readRDS(file="../output/NeiGWAS_Wuest_et_al_sim_all_biomassZ.rds")),
GIF(d=readRDS(file="../output/NeiGWAS_Wuest_et_al_sim_all_floweringZ.rds")),
GIF(d=readRDS(file="../output/NeiGWAS_Wuest_et_al_sim_all_rosetteZ.rds"))
), 3,3)

mgif = as.data.frame(mgif)

colnames(mgif) = c("self","nei","selfxnei")
rownames(mgif) = c("biomass","flowering","rosette")

write.csv(mgif,"./output/Wuest_GIF.csv")

