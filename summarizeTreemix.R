source("plotting_funcs.R")
pve<-matrix(NA,nrow=23,ncol=9)
mig<-seq(0,8)
sds<-matrix(NA,nrow=23,ncol=9)
eigenConc<-matrix(NA,nrow=23,ncol=9)
pdf("covResidDump.pdf",width=6,height=6)
for(j in 1:23){
	for(i in mig){
	    prefix<-paste("tro_ch",j,"_m",i,sep="")
	    R<-plot_resid(prefix,"poplist.txt")
	    W<-plot_cov(prefix,"poplist.txt")
	    rr<-R[upper.tri(R)]
	    sds[j,i+1]<-sd(rr)
	    pco<-prcomp(R,center=TRUE)
	    pp<-pco$sdev^2/sum(pco$sdev^2)
	    eigenConc[j,i+1]<-sum(pp^2)
	    n<-23
	    Rbar<-sum(R[upper.tri(R)])/(n*(n-1)/2)
	    Wbar<-sum(W[upper.tri(W)])/(n*(n-1)/2)
	    num<-sum((R[upper.tri(R)]-Rbar)^2)
	    dnom<-sum((W[upper.tri(W)]-Wbar)^2)
	    pve[j,i+1]<-1-(num/dnom)
	    plot_tree(prefix)

	}
}
dev.off()

eigenRank<-1/eigenConc
## does require Z
cor.test(pve[,1],eigenRank[,1])
#
#	Pearson's product-moment correlation
#
#data:  pve[, 1] and eigenRank[, 1]
#t = 3.894, df = 21, p-value = 0.0008369
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.3210121 0.8364687
#sample estimates:
#      cor
#0.6475317
cor.test(pve[,1],eigenConc[,1])
#
#	Pearson's product-moment correlation
#
#data:  pve[, 1] and eigenConc[, 1]
#t = -3.1404, df = 21, p-value = 0.00494
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.7927705 -0.1996066
#sample estimates:
#       cor 
#-0.5652971 

pc1<-matrix(NA,nrow=23,ncol=23)
for(j in 1:23){
	i<-0
	prefix<-paste("tro_ch",j,"_m",i,sep="")
        R<-plot_resid(prefix,"poplist.txt")
        pco<-prcomp(R,center=TRUE)
	pc1[,j]<-pco$x[,1]
}
pops<-gsub(pattern="\\d",replacement="",x=colnames(R),perl=TRUE)
rownames(pc1)<-pops

## flip signs to maximize similarity
for(j in 2:23){
	if(cor(pc1[,1],pc1[,j]) < 0){
		pc1[,j]<--1*pc1[,j]
	}
}

## supplement figure PC structure
library(corrgram)

cl<-1.3;ca<-1.1
pdf("sfig_chPc1A.pdf",width=5,height=5)
par(mar=c(4.5,4.5,1,1))
corrgram(pc1,upper.panel=panel.cor,lower.panel=panel.pts,pch=20)
dev.off()
pdf("sfig_chPc1B.pdf",width=5,height=5)
par(mar=c(4.5,4.5,1.5,1.5))
plot(pcc$x[,1],pcc$x[,2],type='n',xlab="PC 1 (52%)",ylab="PC 2 (24%)",cex.lab=cl,cex.axis=ca)
text(pcc$x[,1],pcc$x[,2],c(1:22,"Z"))
title(main="(B) PCA",cex.main=cl)
dev.off()

ll<-matrix(NA,nrow=23,ncol=9)
for(j in 1:23){
	for(i in mig){
		lf<-paste("tro_ch",j,"_m",i,".llik",sep="")
		ld<-read.table(lf)
		ll[j,i+1]<-ld[2,7]
	}
}


pdf("treemixLL.pdf",width=8,height=5)
par(mfrow=c(1,2))
for(j in 1:23){
	plot(0:8,ll[j,],type='l',xlab="Migration",ylab="Log likelihood")
	title(main=paste("Chromosome ",j,sep=""))
	plot(0:8,pve[j,],type='l',xlab="Migration",ylab="PVE",ylim=c(0.7,1))
	title(main=paste("Chromosome ",j,sep=""))
}
dev.off()

save(list=ls(),file="tmx.rdat")

image(t(pve),axes=FALSE,xlab="Admixture edges",ylab="Chromosome")
axis(1,at=c(0:8)/8,(1:9)-1)
axis(2,at=c(0:22)/22,(1:23))
box()

pdf("treemixTrees.pdf",width=8,height=9)
par(mfrow=c(3,3))
par(mar=c(1,1,1,1))
for(j in 1:23){
	for(i in mig){
		prefix<-paste("tro_ch",j,"_m",i,sep="")
    		plot_tree(prefix,cex=.7)
	}
}
dev.off()

