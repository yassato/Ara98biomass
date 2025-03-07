library(tidyverse)
library(rNeighborGWAS)
library(gaston)

#####################
# Wuest et al. (2022)
pheno = read.csv("../pheno/Wuest2022PLoSBiolData/competition.csv")

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
g = read.csv("../pheno/call_method_75_TAIR9_250k.csv.gz",header=TRUE,skip=1)
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
saveRDS(geno,file="../output/RegMap250k.rds",compress=TRUE,version=2) # for permutation test

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

coval = model.matrix(~Block+posAB+CommunityType,pheno2)
sum(is.na(pheno2$biomass))

Z = model.matrix(~Pot.No,pheno2)
ZZ = tcrossprod(Z) # add random effects of pot IDs
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

res = parallel::mcmapply(test_marker_i,1:q,mc.cores=8L)
out = cbind(geno$pos,AF,t(res))
colnames(out) = c("chr","pos","AF","Z_s","Z_n","Z_sim","p_s","p_n","p_sim")

saveRDS(out,file="../output/NeiGWAS_Wuest_et_al_sim_all_biomassZ.rds")
saveRDS(pheno2,file="../output/Wuest_pheno2.rds") # for permutation tests
