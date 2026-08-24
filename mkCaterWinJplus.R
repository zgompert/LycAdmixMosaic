library(ape)
library(phytools)

## 10 kb windows BTBxGNPxSIN 
# A = BTB, B = GNP, C = SIN

scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
for(ch in 1:23){
	inf<-paste("winBTBxGNPxSIN/winout",ch,".tsv",sep="")
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

## (BTB,SIN), GNP
mat<-rbind(unlist(scAC),unlist(scAB),unlist(scBC))

apply(mat,1,mean)
# [1] 0.7360776 0.1463881 0.1175343

szs<-unlist(lapply(scAB,length))
ch<-rep(c(1:23),szs)

matBTBxSINxGNP<-mat

## 10 kb windows
# A = BCR, B = GNP, C = SIN

scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
for(ch in 1:23){
	inf<-paste("winBCRxGNPxSIN/winout",ch,".tsv",sep="")
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

## (BCR,SIN), GNP
mat<-rbind(unlist(scAC),unlist(scAB),unlist(scBC))

apply(mat,1,mean)
# 0.737885 0.143361 0.118754

szs<-unlist(lapply(scAB,length))
ch<-rep(c(1:23),szs)

matBCRxSINxGNP<-mat

## 10 kb windows 
# A = HNV, B = GNP, C = SIN

scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
for(ch in 1:23){
	inf<-paste("winHNVxGNPxSIN/winout",ch,".tsv",sep="")
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

## rearrange (HNV,SIN),GNP
mat<-rbind(unlist(scAC),unlist(scAB),unlist(scBC))

apply(mat,1,mean)
#[1] [1] 0.4881710 0.2766805 0.2351485

szs<-unlist(lapply(scAB,length))
ch<-rep(c(1:23),szs)

matHNVxSINxGNP<-mat


## 10 kb windows 
# A = ABM, B = SIN, C = TBY

scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
for(ch in 1:23){
	inf<-paste("winABMxSINxTBY/winout",ch,".tsv",sep="")
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
#[1] 0.81096749 0.09340857 0.09562394

szs<-unlist(lapply(scAB,length))
ch<-rep(c(1:23),szs)

matABMxSINxTBY<-mat

## 10 kb windows 
# A = BHP, B = SIN, C = YG

scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
for(ch in 1:23){
	inf<-paste("winBHPxSINxYG/winout",ch,".tsv",sep="")
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
#[1] ] 0.83663213 0.09973298 0.06363489

szs<-unlist(lapply(scAB,length))
ch<-rep(c(1:23),szs)

matBHPxSINxYG<-mat

## 10 kb windows 
# A = BCR, B = BTB, C = YG 

scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
for(ch in 1:23){
	inf<-paste("winBCRxBTBxYG/winout",ch,".tsv",sep="")
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
#[1]  0.84042423 0.07961246 0.07996331

szs<-unlist(lapply(scAB,length))
ch<-rep(c(1:23),szs)

matBCRxBTBxYG<-mat

## 10 kb windows 
# A = LS, B = YG, C = SIN

scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
for(ch in 1:23){
	inf<-paste("winLSxYGxSIN/winout",ch,".tsv",sep="")
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
#[1] 0.96683269 0.01797436 0.01519296

szs<-unlist(lapply(scAB,length))
ch<-rep(c(1:23),szs)

matLSxYGxSIN<-mat

ccs<-c("darkslategray","lightsalmon3","honeydew3")


cl<-1.4;ca<-1.1;ct<-1.05;cm<-1.3
pdf("fig_caterWinJplus.pdf",width=8,height=10)
layout(matrix(1:14,nrow=7,ncol=2,byrow=TRUE),widths=c(1.5,7.5),heights=rep(2,7))

par(mar=c(1,2.5,1,1))
ttxt<-"((BTB,SIN),GNP);"
tre<-read.tree(text=ttxt)
plot.phylo(tre,type="cladogram",direction="upwards",tip.color=ccs[1],edge.color=ccs[1],node.color=ccs[1],cex=ct,label.offset=.1)
title(main="(A)",cex.main=cm)

par(mar=c(2,4.5,1,1))
ccs<-c("darkslategray","lightsalmon3","honeydew3")
xx<-barplot(matBTBxSINxGNP,border=NA,axes=FALSE,xlab="",ylab="Score",col=ccs,cex.lab=cl,cex.axis=ca)
mids<-tapply(X=xx,INDEX=ch,mean)
bnds<-tapply(X=xx,INDEX=ch,max)[-23]
abline(v=bnds)
axis(2)
axis(1,mids,c(1:22,"Z"))

par(mar=c(1,2.5,1,1))
ttxt<-"((BCR,SIN),GNP);"
tre<-read.tree(text=ttxt)
plot.phylo(tre,type="cladogram",direction="upwards",tip.color=ccs[1],edge.color=ccs[1],node.color=ccs[1],cex=ct,label.offset=.1)
title(main="(B)",cex.main=cm)

par(mar=c(2,4.5,1,1))
ccs<-c("darkslategray","lightsalmon3","honeydew3")
xx<-barplot(matBCRxSINxGNP,border=NA,axes=FALSE,xlab="",ylab="Score",col=ccs,cex.lab=cl,cex.axis=ca)
mids<-tapply(X=xx,INDEX=ch,mean)
bnds<-tapply(X=xx,INDEX=ch,max)[-23]
abline(v=bnds)
axis(2)
axis(1,mids,c(1:22,"Z"))

par(mar=c(1,2.5,1,1))
ttxt<-"((HNV,SIN),GNP);"
tre<-read.tree(text=ttxt)
plot.phylo(tre,type="cladogram",direction="upwards",tip.color=ccs[1],edge.color=ccs[1],node.color=ccs[1],cex=ct,label.offset=.1)
title(main="(C)",cex.main=cm)

par(mar=c(2,4.5,1,1))
ccs<-c("darkslategray","lightsalmon3","honeydew3")
xx<-barplot(matHNVxSINxGNP,border=NA,axes=FALSE,xlab="",ylab="Score",col=ccs,cex.lab=cl,cex.axis=ca)
mids<-tapply(X=xx,INDEX=ch,mean)
bnds<-tapply(X=xx,INDEX=ch,max)[-23]
abline(v=bnds)
axis(2)
axis(1,mids,c(1:22,"Z"))

par(mar=c(1,2.5,1,1))
ttxt<-"((ABM,SIN),TBY);"
tre<-read.tree(text=ttxt)
plot.phylo(tre,type="cladogram",direction="upwards",tip.color=ccs[1],edge.color=ccs[1],node.color=ccs[1],cex=ct,label.offset=.1)
title(main="(D)",cex.main=cm)

par(mar=c(2,4.5,1,1))

ccs<-c("darkslategray","lightsalmon3","honeydew3")
xx<-barplot(matABMxSINxTBY,border=NA,axes=FALSE,xlab="",ylab="Score",col=ccs,cex.lab=cl,cex.axis=ca)
mids<-tapply(X=xx,INDEX=ch,mean)
bnds<-tapply(X=xx,INDEX=ch,max)[-23]
abline(v=bnds)
axis(2)
axis(1,mids,c(1:22,"Z"))

par(mar=c(1,2.5,1,1))
ttxt<-"((BHP,SIN),YG);"
tre<-read.tree(text=ttxt)
plot.phylo(tre,type="cladogram",direction="upwards",tip.color=ccs[1],edge.color=ccs[1],node.color=ccs[1],cex=ct,label.offset=.1)
title(main="(E)",cex.main=cm)

par(mar=c(2,4.5,1,1))

ccs<-c("darkslategray","lightsalmon3","honeydew3")
xx<-barplot(matBHPxSINxYG,border=NA,axes=FALSE,xlab="",ylab="Score",col=ccs,cex.lab=cl,cex.axis=ca)
mids<-tapply(X=xx,INDEX=ch,mean)
bnds<-tapply(X=xx,INDEX=ch,max)[-23]
abline(v=bnds)
axis(2)
axis(1,mids,c(1:22,"Z"))

par(mar=c(1,2.5,1,1))
ttxt<-"((BCR,BTB),YG);"
tre<-read.tree(text=ttxt)
plot.phylo(tre,type="cladogram",direction="upwards",tip.color=ccs[1],edge.color=ccs[1],node.color=ccs[1],cex=ct,label.offset=.1)
title(main="(F)",cex.main=cm)

par(mar=c(2,4.5,1,1))

ccs<-c("darkslategray","lightsalmon3","honeydew3")
xx<-barplot(matBCRxBTBxYG,border=NA,axes=FALSE,xlab="",ylab="Score",col=ccs,cex.lab=cl,cex.axis=ca)
mids<-tapply(X=xx,INDEX=ch,mean)
bnds<-tapply(X=xx,INDEX=ch,max)[-23]
abline(v=bnds)
axis(2)
axis(1,mids,c(1:22,"Z"))

par(mar=c(1,2.5,1,1))
ttxt<-"((LS,YG),SIN);"
tre<-read.tree(text=ttxt)
plot.phylo(tre,type="cladogram",direction="upwards",tip.color=ccs[1],edge.color=ccs[1],node.color=ccs[1],cex=ct,label.offset=.1)
title(main="(G)",cex.main=cm)

par(mar=c(2,4.5,1,1))

ccs<-c("darkslategray","lightsalmon3","honeydew3")
xx<-barplot(matLSxYGxSIN,border=NA,axes=FALSE,xlab="",ylab="Score",col=ccs,cex.lab=cl,cex.axis=ca)
mids<-tapply(X=xx,INDEX=ch,mean)
bnds<-tapply(X=xx,INDEX=ch,max)[-23]
abline(v=bnds)
axis(2)
axis(1,mids,c(1:22,"Z"))

dev.off()
