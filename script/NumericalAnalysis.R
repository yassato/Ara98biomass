# Numerical analysis and figures for FDS and complementarity
library(tidyverse)
library(patchwork)

##########
# dominant case
W_AA = function(x,b2,b12) (b12+b2)*(4*x-2*x^2-1)+1
W_Aa = function(x,b2,b12) (b12+b2)*(4*x-2*x^2-1)+1
W_aa = function(x,b2,b12) (b12-b2)*(4*x-2*x^2-1)+1

W_A = function(x,b2,b12) x*W_AA(x,b2,b12)+(1-x)*W_Aa(x,b2,b12)
W_a = function(x,b2,b12) x*W_Aa(x,b2,b12)+(1-x)*W_aa(x,b2,b12)
W_m = function(x,b2,b12) x*W_A(x,b2,b12)+(1-x)*W_a(x,b2,b12)

dRY_AA = function(x,b2,b12) (x^2)*(W_AA(x,b2,b12)/W_AA(1,b2,b12)) - (x^2)*W_AA(1,b2,b12)
dRY_Aa = function(x,b2,b12) (2*x*(1-x))*(W_Aa(x,b2,b12)/W_Aa(1,b2,b12)) - (2*x*(1-x))*W_Aa(1,b2,b12)
dRY_aa = function(x,b2,b12) ((1-x)^2)*(W_aa(x,b2,b12)/W_aa(0,b2,b12)) - ((1-x)^2)*W_aa(0,b2,b12)

compl = function(x,b2,b12) 3*((dRY_AA(x,b2,b12) + dRY_Aa(x,b2,b12) + dRY_aa(x,b2,b12))/3)*((W_AA(1,b2,b12)+W_Aa(0,b2,b12)+W_aa(0,b2,b12))/3)
sel = function(x,b2,b12) {
  M_bar = (W_AA(1,b2,b12)+W_Aa(1,b2,b12)+W_aa(0,b2,b12))/3
  RY_bar = (dRY_AA(x,b2,b12)+dRY_Aa(x,b2,b12)+dRY_aa(x,b2,b12))/3
  
  covAA = (dRY_AA(x,b2,b12) - RY_bar)*(W_AA(1,b2,b12) - M_bar)
  covAa = (dRY_Aa(x,b2,b12) - RY_bar)*(W_Aa(1,b2,b12) - M_bar)
  covaa = (dRY_aa(x,b2,b12) - RY_bar)*(W_aa(0,b2,b12) - M_bar)
  return(2*(covAA+covAa+covaa))
}
bioe = function(x,b2,b12) compl(x,b2,b12) + sel(x,b2,b12)

wplot = function(b2,b12,pch) {
  p = ggplot(NULL,aes(x=0.293,y=0.293*W_A(0.293,b2,b12)+(1-0.293)*W_a(0.293,b2,b12))) + 
    geom_point(pch=pch,size=3) + ylim(0.5,1.5) + xlim(0,1) + theme_classic() +
    geom_function(fun=W_A, col=grey(0.0,1.0),args=list(b2,b12),lty=2) +
    geom_function(fun=W_a, col=grey(0.0,0.5),args=list(b2,b12),lty=2) + 
    geom_function(fun=W_m, lty=1,args=list(b2,b12)) + 
    ylab("") + xlab("") +
    labs(title=bquote(beta[2]*" = "*.(b2)*"; "*beta[12]*" = "*.(b12)))
  return(p)
} 

pa = wplot(b2=-0.3,b12=0.0,pch=16)
pa = pa + geom_text(aes(x=0.025,y=1.5),label="- - - (black): A alleles",size=3,hjust=0) + 
  geom_text(aes(x=0.025,y=1.4),label="- - - (gray): a alleles",colour=grey(0.5,1),size=3,hjust=0) +
  geom_text(aes(x=0.025,y=1.3),label="-----  (solid): mean",size=3,hjust=0) + 
  ylab("Absolute fitness")
pb = wplot(b2=-0.3,b12=-0.2,pch=16)
pc = wplot(b2=-0.3,b12=0.2,pch=16)
pd = wplot(b2=0.3,b12=0.0,pch=1)
pe = wplot(b2=0.3,b12=-0.2,pch=1)
pf = wplot(b2=0.3,b12=0.2,pch=1)

bplot = function(b2,b12,pch) {
  p = ggplot(NULL,aes(x=0.293,y=compl(0.293,b2,b12))) + 
    geom_point(pch=pch,size=3) + ylim(-1.5,1.5) + xlim(0,1) + theme_classic() +
    geom_function(fun=compl,col=rgb(0,0,1,1),args=list(b2,b12)) +
    geom_function(fun=sel,col=rgb(1,0,0,1),args=list(b2,b12)) + 
    geom_function(fun=bioe,lwd=1.25,col=grey(0.5,0.5),args=list(b2,b12)) + 
    ylab("") + xlab("Frequency of A alleles") +
    labs(title=bquote(beta[2]*" = "*.(b2)*"; "*beta[12]*" = "*.(b12)))
  return(p)
} 

ba = bplot(b2=-0.3,b12=0.0,pch=16)
ba = ba + geom_text(aes(x=0.025,y=1.25),label="----- (blue): complementarity",colour=rgb(0,0,1,1),size=3,hjust=0) + 
  geom_text(aes(x=0.025,y=1.0),label="----- (red): selection",colour=rgb(1,0,0,1),size=3,hjust=0) +
  geom_text(aes(x=0.025,y=0.75),label="----- (grey): net effect",colour=grey(0.25,0.75),size=3,hjust=0) +
  ylab("Biodiversity effects")
bb = bplot(b2=-0.3,b12=-0.2,pch=16)
bc = bplot(b2=-0.3,b12=0.2,pch=16)
bd = bplot(b2=0.3,b12=0.0,pch=1)
be = bplot(b2=0.3,b12=-0.2,pch=1)
bf = bplot(b2=0.3,b12=0.2,pch=1)

bp = ((pa | pb | pc | pd | pe | pf) / (ba | bb | bc | bd | be | bf)) + plot_annotation(tag_levels = "a")
ggsave(bp,filename="./NumericalAnalysis/dominant.pdf",width=16,height=6)

len = length(seq(-0.4,0.4,by=0.01))
X = rep(seq(-0.4,0.4,by=0.01),each=len)
Y = rep(seq(-0.4,0.4,by=0.01),len)
C = mapply(compl,x=rep(0.293,len),b2=X,b12=Y)
S = mapply(sel,x=rep(0.293,len),b2=X,b12=Y)
p_comp = ggplot(NULL,aes(x=Y,y=X,fill=C)) + 
  ylab(expression(beta[2])) + xlab(expression(beta[12])) +
  geom_raster() + theme_classic() + 
  labs(title="Complementarity effects") +
  scale_fill_gradient2(name="",low="navy", mid="white", high="red", midpoint=0)

p_sel = ggplot(NULL,aes(x=Y,y=X,fill=S)) + 
  ylab(expression(beta[2])) + xlab(expression(beta[12])) +
  geom_raster() + theme_classic() + 
  labs(title="Selection effects") +
  scale_fill_gradient2(name="",low="navy", mid="white", high="red", midpoint=0)

pbio = p_comp + p_sel
ggsave(pbio,filename="./NumericalAnalysis/biodominant.pdf",width=9,height=4)


############
# additive case
W_AA = function(x,b2,b12) (b12+b2)*(2*x-1)+1
W_Aa = function(x,b2,b12) 1
W_aa = function(x,b2,b12) (b12-b2)*(2*x-1)+1

dRY_AA = function(x,b2,b12) (x^2)*(W_AA(x,b2,b12)/W_AA(1,b2,b12)) - (x^2)*W_AA(1,b2,b12)
dRY_Aa = function(x,b2,b12) (2*x*(1-x))*(W_Aa(x,b2,b12)/W_Aa(1,b2,b12)) - (2*x*(1-x))*W_Aa(1,b2,b12)
dRY_aa = function(x,b2,b12) ((1-x)^2)*(W_aa(x,b2,b12)/W_aa(0,b2,b12)) - ((1-x)^2)*W_aa(0,b2,b12)

wplot = function(b2,b12,pch=c(1,16)) {
  p2 = 0.5*(b12-b2)/b12
  w2 = p2*W_A(p2,b2,b12)+(1-p2)*W_a(p2,b2,b12)
  p = ggplot(NULL,aes(x=0.5,y=0.5*W_A(0.5,b2,b12)+(1-0.5)*W_a(0.5,b2,b12))) + 
    geom_point(pch=pch[1],size=3) + ylim(0.5,1.5) + xlim(0,1) + theme_classic() +
    geom_function(fun=W_A, col=grey(0.0,1.0),args=list(b2,b12)) +
    geom_function(fun=W_a, col=grey(0.0,0.5),args=list(b2,b12)) + 
    geom_function(fun=W_m, lty=2,args=list(b2,b12)) + 
    ylab("Absolute fitness") + xlab("Frequency of A alleles") +
    labs(title=bquote(beta[2]*" = "*.(b2)*"; "*beta[12]*" = "*.(b12))) + 
    geom_point(aes(x=p2,y=w2),pch=pch[2],size=3)
  return(p)
} 

pa = wplot(b2=-0.3,b12=0,pch=c(16,1))
pa = pa + geom_text(aes(x=0.025,y=1.5),label="----- (black): A alleles",size=3,hjust=0) + 
  geom_text(aes(x=0.025,y=1.4),label="----- (gray): a alleles",colour=grey(0.5,1),size=3,hjust=0) +
  geom_text(aes(x=0.025,y=1.3),label="- - -  (dashed): mean",size=3,hjust=0)
pb = wplot(b2=-0.3,b12=-0.2,pch=c(16,1))
pc = wplot(b2=-0.3,b12=0.2,pch=c(16,1))
pd = wplot(b2=0.3,b12=0,pch=c(1,16))
pe = wplot(b2=0.3,b12=-0.2,pch=c(1,16))
pf = wplot(b2=0.3,b12=0.2,pch=c(1,16))

bplot = function(b2,b12,pch=c(1,16)) {
  p2 = 0.5*(b12-b2)/b12
  p = ggplot(data=data_frame(x=c(0.5,p2),y=c(compl(0.5,b2,b12),compl(p2,b2,b12))),aes(x=x,y=y)) + 
    geom_point(pch=pch,size=3) + ylim(-1.5,1.5) + xlim(0,1) + theme_classic() +
    geom_function(fun=compl,col=rgb(0,0,1,1),args=list(b2,b12)) +
    geom_function(fun=sel,col=rgb(1,0,0,1),args=list(b2,b12)) + 
    geom_function(fun=bioe,lwd=1.25,col=grey(0.5,0.5),args=list(b2,b12)) + 
    ylab("") + xlab("Frequency of A alleles") +
    labs(title=bquote(beta[2]*" = "*.(b2)*"; "*beta[12]*" = "*.(b12)))
  return(p)
} 

ba = bplot(b2=-0.3,b12=0,pch=c(16,1))
ba = ba + geom_text(aes(x=0.025,y=1.25),label="----- (blue): complementarity",colour=rgb(0,0,1,1),size=3,hjust=0) + 
  geom_text(aes(x=0.025,y=1.0),label="----- (red): selection",colour=rgb(1,0,0,1),size=3,hjust=0) +
  geom_text(aes(x=0.025,y=0.75),label="----- (grey): net effect",colour=grey(0.25,0.75),size=3,hjust=0) +
  ylab("Biodiversity effects")
bb = bplot(b2=-0.3,b12=-0.2,pch=c(16,1))
bc = bplot(b2=-0.3,b12=0.2,pch=c(16,1))
bd = bplot(b2=0.3,b12=0,pch=c(1,16))
be = bplot(b2=0.3,b12=-0.2,pch=c(1,16))
bf = bplot(b2=0.3,b12=0.2,pch=c(1,16))

bp = ((pa | pb | pc | pd | pe | pf) / (ba | bb | bc | bd | be | bf)) + plot_annotation(tag_levels = "a")
ggsave(bp,filename="./NumericalAnalysis/additive.pdf",width=16,height=6)

len = length(seq(-0.4,0.4,by=0.01))
X = rep(seq(-0.4,0.4,by=0.01),each=len)
Y = rep(seq(-0.4,0.4,by=0.01),len)
C = mapply(compl,x=rep(0.5,len),b2=X,b12=Y)
S = mapply(sel,x=rep(0.5,len),b2=X,b12=Y)
p_comp = ggplot(NULL,aes(x=Y,y=X,fill=C)) + 
  ylab(expression(beta[2])) + xlab(expression(beta[12])) +
  geom_raster() + theme_classic() + 
  labs(title="Complementarity effects") +
  scale_fill_gradient2(name="",low="navy", mid="white", high="red", midpoint=0)

p_sel = ggplot(NULL,aes(x=Y,y=X,fill=S)) + 
  ylab(expression(beta[2])) + xlab(expression(beta[12])) +
  geom_raster() + theme_classic() + 
  labs(title="Selection effects") +
  scale_fill_gradient2(name="",low="navy", mid="white", high="red", midpoint=0)

pbio = p_comp + p_sel
ggsave(pbio,filename="./NumericalAnalysis/bioadditive.pdf",width=9,height=4)


###########
# inbred case
W_A = function(x,b2,b12) (b12+b2)*(2*x-1)+1
W_a = function(x,b2,b12) (b12-b2)*(2*x-1)+1
W_m = function(x,b2,b12) 2*b2*x*(2*x-1)+(b12-b2)*(2*x-1)+1

dRY_A = function(x,b2,b12) x*(W_A(x,b2,b12)/W_A(1,b2,b12)) - x*W_A(1,b2,b12)
dRY_a = function(x,b2,b12) (1-x)*(W_a(x,b2,b12)/W_a(0,b2,b12)) - (1-x)*W_a(0,b2,b12)

compl = function(x,b2,b12) 2*((dRY_A(x,b2,b12) + dRY_a(x,b2,b12))/2)*((W_A(1,b2,b12)+W_a(0,b2,b12))/2)
sel = function(x,b2,b12) {
  M_bar = (W_A(1,b2,b12)+W_a(0,b2,b12))/2
  RY_bar = (dRY_A(x,b2,b12)+dRY_a(x,b2,b12))/2
  
  covA = (dRY_A(x,b2,b12) - RY_bar)*(W_A(1,b2,b12) - M_bar)
  cova = (dRY_a(x,b2,b12) - RY_bar)*(W_a(0,b2,b12) - M_bar)
  return(2*(covA+cova))
}
bioe = function(x,b2,b12) compl(x,b2,b12) + sel(x,b2,b12)


wplot = function(b2,b12,pch) {
  p = ggplot(NULL,aes(x=0.5,y=1)) + 
    geom_point(pch=pch,size=3) + ylim(0.5,1.5) + xlim(0,1) + theme_classic() +
    geom_function(fun=function(x,b2,b12) (b12+b2)*(2*x-1)+1, col=grey(0.0,1.0),args=list(b2,b12)) +
    geom_function(fun=function(x,b2,b12) (b12-b2)*(2*x-1)+1, col=grey(0.0,0.5),args=list(b2,b12)) + 
    geom_function(fun=function(x,b2,b12) 2*b2*x*(2*x-1)+(b12-b2)*(2*x-1)+1, lty=2,args=list(b2,b12)) + 
    ylab("") + xlab("") +
    labs(title=bquote(beta[2]*" = "*.(b2)*"; "*beta[12]*" = "*.(b12)))
  return(p)
} 

pa = wplot(b2=-0.3,b12=0.0,pch=16)
pa = pa + geom_text(aes(x=0.025,y=1.5),label="----- (black): A alleles",size=3,hjust=0) + 
  geom_text(aes(x=0.025,y=1.4),label="----- (gray): a alleles",colour=grey(0.5,1),size=3,hjust=0) +
  geom_text(aes(x=0.025,y=1.3),label="- - -  (dashed): mean",size=3,hjust=0) +
  ylab("Absolute fitness")
pb = wplot(b2=-0.3,b12=-0.2,pch=16)
pc = wplot(b2=-0.3,b12=0.2,pch=16)
pd = wplot(b2=0.3,b12=0.0,pch=1)
pe = wplot(b2=0.3,b12=-0.2,pch=1)
pf = wplot(b2=0.3,b12=0.2,pch=1)

bplot = function(b2,b12,pch) {
  p = ggplot(NULL,aes(x=0.5,y=compl(0.5,b2,b12))) + 
    geom_point(pch=pch,size=3) + ylim(-1.5,1.5) + xlim(0,1) + theme_classic() +
    geom_function(fun=compl,col=rgb(0,0,1,1),args=list(b2,b12)) +
    geom_function(fun=sel,col=rgb(1,0,0,1),args=list(b2,b12)) + 
    geom_function(fun=bioe,lwd=1.25,col=grey(0.5,0.5),args=list(b2,b12)) + 
    ylab("") + xlab("Frequency of A alleles") +
    labs(title=bquote(beta[2]*" = "*.(b2)*"; "*beta[12]*" = "*.(b12)))
  return(p)
} 

ba = bplot(b2=-0.3,b12=0.0,pch=16)
ba = ba + geom_text(aes(x=0.025,y=1.25),label="----- (blue): complementarity",colour=rgb(0,0,1,1),size=3,hjust=0) + 
  geom_text(aes(x=0.025,y=1.0),label="----- (red): selection",colour=rgb(1,0,0,1),size=3,hjust=0) +
  geom_text(aes(x=0.025,y=0.75),label="----- (grey): net effect",colour=grey(0.25,0.75),size=3,hjust=0) +
  ylab("Biodiversity effects")
bb = bplot(b2=-0.3,b12=-0.2,pch=16)
bc = bplot(b2=-0.3,b12=0.2,pch=16)
bd = bplot(b2=0.3,b12=0.0,pch=1)
be = bplot(b2=0.3,b12=-0.2,pch=1)
bf = bplot(b2=0.3,b12=0.2,pch=1)

bp = ((pa | pb | pc | pd | pe | pf) / (ba | bb | bc | bd | be | bf)) + plot_annotation(tag_levels = "a")
ggsave(bp,filename="./NumericalAnalysis/inbred.pdf",width=16,height=6)

len = length(seq(-0.4,0.4,by=0.01))
X = rep(seq(-0.4,0.4,by=0.01),each=len)
Y = rep(seq(-0.4,0.4,by=0.01),len)
C = mapply(compl,x=rep(0.5,len),b2=X,b12=Y)
S = mapply(sel,x=rep(0.5,len),b2=X,b12=Y)
p_comp = ggplot(NULL,aes(x=Y,y=X,fill=C)) + 
  ylab(expression(beta[2])) + xlab(expression(beta[12])) +
  geom_raster() + theme_classic() + 
  labs(title="Complementarity effects") +
  scale_fill_gradient2(name="",low="navy", mid="white", high="red", midpoint=0)

p_sel = ggplot(NULL,aes(x=Y,y=X,fill=S)) + 
  ylab(expression(beta[2])) + xlab(expression(beta[12])) +
  geom_raster() + theme_classic() + 
  labs(title="Selection effects") +
  scale_fill_gradient2(name="",low="navy", mid="white", high="red", midpoint=0)

pbio = p_comp + p_sel
ggsave(pbio,filename="./NumericalAnalysis/bioinbred.pdf",width=9,height=4)

