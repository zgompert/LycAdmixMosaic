library(ape)
library(phytools)

## mean within chromsome variance
mnvar<-rep(NA,7)
## var in chromosomem eans
varmn<-rep(NA,7)

## 10 kb windows CLHxTICxBHP
# A = CLH, B = TIC, C = BHP

scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
for(ch in 1:23){
	inf<-paste("winCLHxTICxBHP/winout",ch,".tsv",sep="")
	dat<-read.table(inf,header=TRUE)
	L<-dim(dat)[1]
	NL<-L-5
	ab<-rep(NA,NL)
	ac<-rep(NA,NL)
	bc<-rep(NA,NL)
	for(i in 1:NL){
		ab[i]<-mean(dat$A.B[i:(i+4)])
		bc[i]<-mean(dat$B.C[i:(i+4)])
		ac[i]<-mean(dat$A.C[i:(i+4)])
		sc<-c(ab[i],bc[i],ac[i])
		ssc<-sum(sc)
		ab[i]<-ab[i]/ssc
		bc[i]<-bc[i]/ssc
		ac[i]<-ac[i]/ssc
	}
	scAB[[ch]]<-ab
	scAC[[ch]]<-ac
	scBC[[ch]]<-bc
}

mat<-rbind(unlist(scAB),unlist(scAC),unlist(scBC))

apply(mat,1,mean)
#[1] 0.3816334 0.4418192 0.1765475
# more of genome goes CLH x BHP than CLH x TIC, but Z goes strongly CLH x TIC

szs<-unlist(lapply(scAB,length))
ch<-rep(c(1:23),szs)

mnvar[1]<-(mean(tapply(X=mat[1,],INDEX=ch,var))+mean(tapply(X=mat[2,],INDEX=ch,var))+mean(tapply(X=mat[3,],INDEX=ch,var)))/3
varmn[1]<-(var(tapply(X=mat[1,],INDEX=ch,mean))+var(tapply(X=mat[2,],INDEX=ch,mean))+var(tapply(X=mat[3,],INDEX=ch,mean)))/3
matCLHxTICxBHP<-mat

## 10 kb windows
# A = CP, B = YB, C = TIC

scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
for(ch in 1:23){
	inf<-paste("winCPxYGxTIC/winout",ch,".tsv",sep="")
	dat<-read.table(inf,header=TRUE)
	L<-dim(dat)[1]
	NL<-L-5
	ab<-rep(NA,NL)
	ac<-rep(NA,NL)
	bc<-rep(NA,NL)
	for(i in 1:NL){
		ab[i]<-mean(dat$A.B[i:(i+4)])
		bc[i]<-mean(dat$B.C[i:(i+4)])
		ac[i]<-mean(dat$A.C[i:(i+4)])
		sc<-c(ab[i],bc[i],ac[i])
		ssc<-sum(sc)
		ab[i]<-ab[i]/ssc
		bc[i]<-bc[i]/ssc
		ac[i]<-ac[i]/ssc
	}
	scAB[[ch]]<-ab
	scAC[[ch]]<-ac
	scBC[[ch]]<-bc
}

## rearrange CP,TIC YG
mat<-rbind(unlist(scAC),unlist(scAB),unlist(scBC))

apply(mat,1,mean)
#[1] 0.67438798 0.25421605 0.07139598

szs<-unlist(lapply(scAB,length))
ch<-rep(c(1:23),szs)

mnvar[2]<-(mean(tapply(X=mat[1,],INDEX=ch,var))+mean(tapply(X=mat[2,],INDEX=ch,var))+mean(tapply(X=mat[3,],INDEX=ch,var)))/3
varmn[2]<-(var(tapply(X=mat[1,],INDEX=ch,mean))+var(tapply(X=mat[2,],INDEX=ch,mean))+var(tapply(X=mat[3,],INDEX=ch,mean)))/3

mns<-tapply(X=mat[1,],INDEX=ch,mean)
cor.test(mns[-23],szs[-23])
#t = 8.726, df = 20, p-value = 2.97e-08
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.7495596 0.9537066
#sample estimates:
#      cor
#0.8899313
matCPxTICxYG<-mat



## 10 kb windows 
# A = MR, B = CP, C = YG

scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
for(ch in 1:23){
	inf<-paste("winMRxCPxYG/winout",ch,".tsv",sep="")
	dat<-read.table(inf,header=TRUE)
	L<-dim(dat)[1]
	NL<-L-5
	ab<-rep(NA,NL)
	ac<-rep(NA,NL)
	bc<-rep(NA,NL)
	for(i in 1:NL){
		ab[i]<-mean(dat$A.B[i:(i+4)])
		bc[i]<-mean(dat$B.C[i:(i+4)])
		ac[i]<-mean(dat$A.C[i:(i+4)])
		sc<-c(ab[i],bc[i],ac[i])
		ssc<-sum(sc)
		ab[i]<-ab[i]/ssc
		bc[i]<-bc[i]/ssc
		ac[i]<-ac[i]/ssc
	}
	scAB[[ch]]<-ab
	scAC[[ch]]<-ac
	scBC[[ch]]<-bc
}

## rearrange MR,YG CP
mat<-rbind(unlist(scAC),unlist(scAB),unlist(scBC))

apply(mat,1,mean)
#[1] 0.5371065 0.3367410 0.1261525

szs<-unlist(lapply(scAB,length))
ch<-rep(c(1:23),szs)

mnvar[3]<-(mean(tapply(X=mat[1,],INDEX=ch,var))+mean(tapply(X=mat[2,],INDEX=ch,var))+mean(tapply(X=mat[3,],INDEX=ch,var)))/3
varmn[3]<-(var(tapply(X=mat[1,],INDEX=ch,mean))+var(tapply(X=mat[2,],INDEX=ch,mean))+var(tapply(X=mat[3,],INDEX=ch,mean)))/3

matMRxYGxCP<-mat
mns<-tapply(X=mat[2,],INDEX=ch,mean)
cor.test(mns[-23],szs[-23])
#t = 5.4502, df = 20, p-value = 2.464e-05
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.5213915 0.9010064
#sample estimates:
#      cor
#0.7730629


## 10 kb windows 
# A = CLH, B = YG, C = BHP

scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
for(ch in 1:23){
	inf<-paste("winCLHxYGxBHP/winout",ch,".tsv",sep="")
	dat<-read.table(inf,header=TRUE)
	L<-dim(dat)[1]
	NL<-L-5
	ab<-rep(NA,NL)
	ac<-rep(NA,NL)
	bc<-rep(NA,NL)
	for(i in 1:NL){
		ab[i]<-mean(dat$A.B[i:(i+4)])
		bc[i]<-mean(dat$B.C[i:(i+4)])
		ac[i]<-mean(dat$A.C[i:(i+4)])
		sc<-c(ab[i],bc[i],ac[i])
		ssc<-sum(sc)
		ab[i]<-ab[i]/ssc
		bc[i]<-bc[i]/ssc
		ac[i]<-ac[i]/ssc
	}
	scAB[[ch]]<-ab
	scAC[[ch]]<-ac
	scBC[[ch]]<-bc
}

mat<-rbind(unlist(scAB),unlist(scAC),unlist(scBC))

apply(mat,1,mean)
#[1] 0.16750530 0.75015411 0.08234059

szs<-unlist(lapply(scAB,length))
ch<-rep(c(1:23),szs)

mnvar[4]<-(mean(tapply(X=mat[1,],INDEX=ch,var))+mean(tapply(X=mat[2,],INDEX=ch,var))+mean(tapply(X=mat[3,],INDEX=ch,var)))/3
varmn[4]<-(var(tapply(X=mat[1,],INDEX=ch,mean))+var(tapply(X=mat[2,],INDEX=ch,mean))+var(tapply(X=mat[3,],INDEX=ch,mean)))/3

matCLHxYGxBHP<-mat

## 10 kb windows 
# A = TCI, B = YG, C = BHP

scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
for(ch in 1:23){
	inf<-paste("winTICxYGxBHP/winout",ch,".tsv",sep="")
	dat<-read.table(inf,header=TRUE)
	L<-dim(dat)[1]
	NL<-L-5
	ab<-rep(NA,NL)
	ac<-rep(NA,NL)
	bc<-rep(NA,NL)
	for(i in 1:NL){
		ab[i]<-mean(dat$A.B[i:(i+4)])
		bc[i]<-mean(dat$B.C[i:(i+4)])
		ac[i]<-mean(dat$A.C[i:(i+4)])
		sc<-c(ab[i],bc[i],ac[i])
		ssc<-sum(sc)
		ab[i]<-ab[i]/ssc
		bc[i]<-bc[i]/ssc
		ac[i]<-ac[i]/ssc
	}
	scAB[[ch]]<-ab
	scAC[[ch]]<-ac
	scBC[[ch]]<-bc
}

mat<-rbind(unlist(scAB),unlist(scAC),unlist(scBC))

apply(mat,1,mean)
#[1] 0.23979493 0.66141942 0.09878566

szs<-unlist(lapply(scAB,length))
ch<-rep(c(1:23),szs)

mnvar[5]<-(mean(tapply(X=mat[1,],INDEX=ch,var))+mean(tapply(X=mat[2,],INDEX=ch,var))+mean(tapply(X=mat[3,],INDEX=ch,var)))/3
varmn[5]<-(var(tapply(X=mat[1,],INDEX=ch,mean))+var(tapply(X=mat[2,],INDEX=ch,mean))+var(tapply(X=mat[3,],INDEX=ch,mean)))/3

matTICxYGxBHP<-mat

## 10 kb windows 
# A = CP, B = YG, C = BHP

scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
for(ch in 1:23){
	inf<-paste("winCPxYGxBHP/winout",ch,".tsv",sep="")
	dat<-read.table(inf,header=TRUE)
	L<-dim(dat)[1]
	NL<-L-5
	ab<-rep(NA,NL)
	ac<-rep(NA,NL)
	bc<-rep(NA,NL)
	for(i in 1:NL){
		ab[i]<-mean(dat$A.B[i:(i+4)])
		bc[i]<-mean(dat$B.C[i:(i+4)])
		ac[i]<-mean(dat$A.C[i:(i+4)])
		sc<-c(ab[i],bc[i],ac[i])
		ssc<-sum(sc)
		ab[i]<-ab[i]/ssc
		bc[i]<-bc[i]/ssc
		ac[i]<-ac[i]/ssc
	}
	scAB[[ch]]<-ab
	scAC[[ch]]<-ac
	scBC[[ch]]<-bc
}

mat<-rbind(unlist(scAB),unlist(scAC),unlist(scBC))

apply(mat,1,mean)
#[1] 0.4142357 0.4908772 0.0948871

szs<-unlist(lapply(scAB,length))
ch<-rep(c(1:23),szs)

mnvar[6]<-(mean(tapply(X=mat[1,],INDEX=ch,var))+mean(tapply(X=mat[2,],INDEX=ch,var))+mean(tapply(X=mat[3,],INDEX=ch,var)))/3
varmn[6]<-(var(tapply(X=mat[1,],INDEX=ch,mean))+var(tapply(X=mat[2,],INDEX=ch,mean))+var(tapply(X=mat[3,],INDEX=ch,mean)))/3

mns<-tapply(X=mat[2,],INDEX=ch,mean)
cor.test(mns[-23],szs[-23])
#t = 8.5014, df = 20, p-value = 4.501e-08
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.7392392 0.9515660
#sample estimates:
#      cor
#0.8850154

matCPxYGxBHP<-mat

## 10 kb windows 
# A = EP, B = SHC, C = BHP

scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
for(ch in 1:23){
	inf<-paste("winEPxSHCxBHP/winout",ch,".tsv",sep="")
	dat<-read.table(inf,header=TRUE)
	L<-dim(dat)[1]
	NL<-L-5
	ab<-rep(NA,NL)
	ac<-rep(NA,NL)
	bc<-rep(NA,NL)
	for(i in 1:NL){
		ab[i]<-mean(dat$A.B[i:(i+4)])
		bc[i]<-mean(dat$B.C[i:(i+4)])
		ac[i]<-mean(dat$A.C[i:(i+4)])
		sc<-c(ab[i],bc[i],ac[i])
		ssc<-sum(sc)
		ab[i]<-ab[i]/ssc
		bc[i]<-bc[i]/ssc
		ac[i]<-ac[i]/ssc
	}
	scAB[[ch]]<-ab
	scAC[[ch]]<-ac
	scBC[[ch]]<-bc
}

mat<-rbind(unlist(scAB),unlist(scAC),unlist(scBC))

apply(mat,1,mean)
#[1]  0.3718465 0.5181919 0.1099616

szs<-unlist(lapply(scAB,length))
ch<-rep(c(1:23),szs)

mnvar[7]<-(mean(tapply(X=mat[1,],INDEX=ch,var))+mean(tapply(X=mat[2,],INDEX=ch,var))+mean(tapply(X=mat[3,],INDEX=ch,var)))/3
varmn[7]<-(var(tapply(X=mat[1,],INDEX=ch,mean))+var(tapply(X=mat[2,],INDEX=ch,mean))+var(tapply(X=mat[3,],INDEX=ch,mean)))/3

matEPxSHCxBHP<-mat

summary(mnvar)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.002765 0.003641 0.005545 0.005452 0.006362 0.009845
summary(varmn)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.001945 0.002641 0.004425 0.004489 0.005750 0.008270

ccs<-c("darkslategray","lightsalmon3","honeydew3")

pdf("fig_exampleTrees.pdf",width=4,height=1.5)
par(mfrow=c(1,3))
par(mar=c(1,2.5,1,1))
ttxt<-"((A,B),C);"
tre<-read.tree(text=ttxt)
plot.phylo(tre,type="cladogram",direction="upwards",tip.color=ccs[1],edge.color=ccs[1],node.color=ccs[1],cex=ct,label.offset=.1)

ttxt<-"((A,C),B);"
tre<-read.tree(text=ttxt)
plot.phylo(tre,type="cladogram",direction="upwards",tip.color=ccs[2],edge.color=ccs[2],node.color=ccs[2],cex=ct,label.offset=.1)

ttxt<-"((B,C),A);"
tre<-read.tree(text=ttxt)
plot.phylo(tre,type="cladogram",direction="upwards",tip.color=ccs[3],edge.color=ccs[3],node.color=ccs[3],cex=ct,label.offset=.1)
dev.off()

cl<-1.4;ca<-1.1;ct<-1.05;cm<-1.3
pdf("fig_caterWinSierra.pdf",width=8,height=10)
layout(matrix(1:14,nrow=7,ncol=2,byrow=TRUE),widths=c(1.5,7.5),heights=rep(2,7))

par(mar=c(1,2.5,1,1))
ttxt<-"((CLH,TIC),BHP);"
tre<-read.tree(text=ttxt)
plot.phylo(tre,type="cladogram",direction="upwards",tip.color=ccs[1],edge.color=ccs[1],node.color=ccs[1],cex=ct,label.offset=.1)
title(main="(A)",cex.main=cm)

par(mar=c(2,4.5,1,1))
ccs<-c("darkslategray","lightsalmon3","honeydew3")
xx<-barplot(matCLHxTICxBHP,border=NA,axes=FALSE,xlab="",ylab="Score",col=ccs,cex.lab=cl,cex.axis=ca)
mids<-tapply(X=xx,INDEX=ch,mean)
bnds<-tapply(X=xx,INDEX=ch,max)[-23]
abline(v=bnds)
axis(2)
axis(1,mids,c(1:22,"Z"))

par(mar=c(1,2.5,1,1))
ttxt<-"((CLH,YG),BHP);"
tre<-read.tree(text=ttxt)
plot.phylo(tre,type="cladogram",direction="upwards",tip.color=ccs[1],edge.color=ccs[1],node.color=ccs[1],cex=ct,label.offset=.1)
title(main="(B)",cex.main=cm)

par(mar=c(2,4.5,1,1))
ccs<-c("darkslategray","lightsalmon3","honeydew3")
xx<-barplot(matCLHxYGxBHP,border=NA,axes=FALSE,xlab="",ylab="Score",col=ccs,cex.lab=cl,cex.axis=ca)
mids<-tapply(X=xx,INDEX=ch,mean)
bnds<-tapply(X=xx,INDEX=ch,max)[-23]
abline(v=bnds)
axis(2)
axis(1,mids,c(1:22,"Z"))

par(mar=c(1,2.5,1,1))
ttxt<-"((TIC,YG),BHP);"
tre<-read.tree(text=ttxt)
plot.phylo(tre,type="cladogram",direction="upwards",tip.color=ccs[1],edge.color=ccs[1],node.color=ccs[1],cex=ct,label.offset=.1)
title(main="(C)",cex.main=cm)

par(mar=c(2,4.5,1,1))
ccs<-c("darkslategray","lightsalmon3","honeydew3")
xx<-barplot(matTICxYGxBHP,border=NA,axes=FALSE,xlab="",ylab="Score",col=ccs,cex.lab=cl,cex.axis=ca)
mids<-tapply(X=xx,INDEX=ch,mean)
bnds<-tapply(X=xx,INDEX=ch,max)[-23]
abline(v=bnds)
axis(2)
axis(1,mids,c(1:22,"Z"))

par(mar=c(1,2.5,1,1))
ttxt<-"((CP,TIC),YG);"
tre<-read.tree(text=ttxt)
plot.phylo(tre,type="cladogram",direction="upwards",tip.color=ccs[1],edge.color=ccs[1],node.color=ccs[1],cex=ct,label.offset=.1)
title(main="(D)",cex.main=cm)

par(mar=c(2,4.5,1,1))

ccs<-c("darkslategray","lightsalmon3","honeydew3")
xx<-barplot(matCPxTICxYG,border=NA,axes=FALSE,xlab="",ylab="Score",col=ccs,cex.lab=cl,cex.axis=ca)
mids<-tapply(X=xx,INDEX=ch,mean)
bnds<-tapply(X=xx,INDEX=ch,max)[-23]
abline(v=bnds)
axis(2)
axis(1,mids,c(1:22,"Z"))

par(mar=c(1,2.5,1,1))
ttxt<-"((MR,YG),CP_);"
tre<-read.tree(text=ttxt)
plot.phylo(tre,type="cladogram",direction="upwards",tip.color=ccs[1],edge.color=ccs[1],node.color=ccs[1],cex=ct,label.offset=.1)
title(main="(E)",cex.main=cm)

par(mar=c(2,4.5,1,1))

ccs<-c("darkslategray","lightsalmon3","honeydew3")
xx<-barplot(matMRxYGxCP,border=NA,axes=FALSE,xlab="",ylab="Score",col=ccs,cex.lab=cl,cex.axis=ca)
mids<-tapply(X=xx,INDEX=ch,mean)
bnds<-tapply(X=xx,INDEX=ch,max)[-23]
abline(v=bnds)
axis(2)
axis(1,mids,c(1:22,"Z"))

par(mar=c(1,2.5,1,1))
ttxt<-"((CP,YG),BHP);"
tre<-read.tree(text=ttxt)
plot.phylo(tre,type="cladogram",direction="upwards",tip.color=ccs[1],edge.color=ccs[1],node.color=ccs[1],cex=ct,label.offset=.1)
title(main="(F)",cex.main=cm)

par(mar=c(2,4.5,1,1))

ccs<-c("darkslategray","lightsalmon3","honeydew3")
xx<-barplot(matCPxYGxBHP,border=NA,axes=FALSE,xlab="",ylab="Score",col=ccs,cex.lab=cl,cex.axis=ca)
mids<-tapply(X=xx,INDEX=ch,mean)
bnds<-tapply(X=xx,INDEX=ch,max)[-23]
abline(v=bnds)
axis(2)
axis(1,mids,c(1:22,"Z"))

par(mar=c(1,2.5,1,1))
ttxt<-"((EP,SHC),BHP);"
tre<-read.tree(text=ttxt)
plot.phylo(tre,type="cladogram",direction="upwards",tip.color=ccs[1],edge.color=ccs[1],node.color=ccs[1],cex=ct,label.offset=.1)
title(main="(G)",cex.main=cm)

par(mar=c(2,4.5,1,1))

ccs<-c("darkslategray","lightsalmon3","honeydew3")
xx<-barplot(matEPxSHCxBHP,border=NA,axes=FALSE,xlab="",ylab="Score",col=ccs,cex.lab=cl,cex.axis=ca)
mids<-tapply(X=xx,INDEX=ch,mean)
bnds<-tapply(X=xx,INDEX=ch,max)[-23]
abline(v=bnds)
axis(2)
axis(1,mids,c(1:22,"Z"))

dev.off()
