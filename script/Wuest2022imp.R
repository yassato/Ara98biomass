library(tidyverse)
library(patchwork)
library(rNeighborGWAS)
library(gaston)
source("coord.R")


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

#####################
# Wuest et al. (2022)
pheno = read.csv("./Atha/Wuest2022PLoSBiolData/competition.csv")

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
g = readRDS("../data/SEW2022_sub_snpsMAF5.rds")
pos = readRDS("../data/SEW_sub_posMAF5.rds")
colnames(g) = paste0("X",colnames(g))
g[g==0] = -1 # convert 0,1 alleles into -1,+1

geno = list(pos=pos,g=g)

rm(g); rm(pos)
gc();gc()

geno$g = geno$g[,levels(factor(c(pheno2$focal,pheno2$neighbor)))]

# MAF cutoff 5% again for full set
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

coval = model.matrix(~Block+posAB+CommunityType,pheno2)
sum(is.na(pheno2$biomass))

Z = model.matrix(~Pot.No,pheno2)
ZZ = tcrossprod(Z) # add random effects of pot IDs... not yet fixed
ZZ = as.matrix(Matrix::nearPD(ZZ,maxit=10^6)$mat)

q = ncol(g_self)
K_self = tcrossprod(g_self)
K_self = ((q - 1)/2 + K_self/2)/(q - 1)
K_self = as.matrix(Matrix::nearPD(K_self,maxit=10^6)$mat)
K_nei = tcrossprod(g_nei)/(q - 1)
K_nei = as.matrix(Matrix::nearPD(K_nei,maxit=10^6)$mat)

K_sim <- tcrossprod(g_sim)/(q - 1)
K_sim <- as.matrix(Matrix::nearPD(K_sim, maxit = 10^6)$mat)

Y = scale(pheno2$biomass)
X = coval

aireml <- gaston::lmm.aireml(Y = Y, X = X, K = list(ZZ, K_self), verbose = FALSE)
Ks <- aireml$tau[1] * ZZ + aireml$tau[2] * K_self
eiKs <- eigen(K_self)
res00 <- gaston::lmm.diago(Y = Y, X = X, eigenK = eiKs, verbose = FALSE)
LL00 <- gaston::lmm.diago.profile.likelihood(tau = res00$tau, s2 = res00$sigma2, Y = Y, X = X, eigenK = eiKs)[1,1]

aireml <- gaston::lmm.aireml(Y = Y, X = X, K = list(ZZ, K_self, K_nei), verbose = FALSE)
K <- aireml$tau[1] * ZZ + aireml$tau[2] * K_self + aireml$tau[3] * K_nei
eiK <- eigen(K)
res01 <- gaston::lmm.diago(Y = Y, X = X, eigenK = eiK, verbose = FALSE)

aireml <- gaston::lmm.aireml(Y = Y, X = X, K = list(ZZ, K_self, K_nei, K_sim), verbose = FALSE)
Ksn <- aireml$tau[1] * ZZ + aireml$tau[2] * K_self + aireml$tau[3] * K_nei + aireml$tau[4] * K_sim
eiKsn <- eigen(Ksn)
res02 <- gaston::lmm.diago(Y = Y, X = X, eigenK = eiKsn, verbose = FALSE)

test_marker_i = function(i) { print(i)
    X0 <- cbind(coval, g_self[, i])
    X1 <- cbind(X0, g_nei[, i])
    X2 <- cbind(X1, g_sim[, i])

    LL_self0 <- gaston::lmm.diago.profile.likelihood(tau = res00$tau, s2 = res00$sigma2, Y = Y, X = X0, eigenK = eiKs)[1,1]
    p_self <- stats::pchisq(-2 * (LL00 - LL_self0), 1, lower.tail = FALSE)
    LL_self1 <- gaston::lmm.diago.profile.likelihood(tau = res01$tau, s2 = res01$sigma2, Y = Y, X = X0, eigenK = eiK)[1,1]
    LL_nei <- gaston::lmm.diago.profile.likelihood(tau = res01$tau, s2 = res01$sigma2, Y = Y, X = X1, eigenK = eiK)[1,1]
    p_nei <- stats::pchisq(-2 * (LL_self1 - LL_nei), 1, lower.tail = FALSE)
    
    LL_nei1 <- gaston::lmm.diago.profile.likelihood(tau = res02$tau, s2 = res02$sigma2, Y = Y, X = X1, eigenK = eiKsn)[1,1]
    LL_sim <- gaston::lmm.diago.profile.likelihood(tau = res02$tau, s2 = res02$sigma2, Y = Y, X = X2, eigenK = eiKsn)[1,1]
    p_sim <- stats::pchisq(-2 * (LL_nei1 - LL_sim), 1, lower.tail = FALSE)
    lmm_nei <- gaston::lmm.diago(Y = Y, X = X2, eigenK = eiKsn, verbose = FALSE)
    z_eff <- lmm_nei$BLUP_beta[-2:0 + length(lmm_nei$BLUP_beta)]/sqrt(diag(lmm_nei$varbeta)[-2:0 + length(lmm_nei$BLUP_beta)])
    resList <- c(z_eff, p_self, p_nei, p_sim)
    return(resList)
}

res = parallel::mcmapply(test_marker_i,1:q,mc.cores=10L)
out = cbind(geno$pos,AF,t(res))
colnames(out) = c("chr","pos","AF","Z_s","Z_n","Z_sim","p_s","p_n","p_sim")

saveRDS(out,file="../output/NeiGWAS_Wuest_et_al_sim_all_biomassZ.rds")


