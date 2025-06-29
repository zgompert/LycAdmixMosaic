## 10 kb windows HJxVExSIN
# A = HJ, B = VE, C = SIN
# intraspecific only? 

scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
for(ch in 1:23){
	inf<-paste("winout",ch,".tsv",sep="")
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
#[1] 0.6354200 0.2473369 0.1172431
# looks like HJ x SIN more than VE x SIN? could be real?


szs<-unlist(lapply(scAB,length))
ch<-rep(c(1:23),szs)

pdf("winHJxVExSIN.pdf",width=9,height=4)
par(mar=c(5,5,1,1))
xx<-barplot(mat,border=NA,axes=FALSE,xlab="Chromosome",ylab="Score")
mids<-tapply(X=xx,INDEX=ch,mean)
bnds<-tapply(X=xx,INDEX=ch,max)[-23]
abline(v=bnds)
axis(2)
axis(1,mids,c(1:22,"Z"))
dev.off()
