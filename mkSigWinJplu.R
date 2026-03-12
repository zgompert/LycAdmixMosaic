## statistical summaries and null model testing for Caster window analyses
## for JH and friends

summarizeWins<-function(M=mat){
	sigs<-apply(M[,3:5],1,sum)
	## defined from tree
	difs<-M[,4]-M[,5]
	## defined from proportions
	os<-order(apply(M[,3:5],2,sum))
	mdifs<-M[,2+os[2]]-M[,2+os[1]]
	## defined from propotion A and Z different
	mnA<-apply(M[M[,1]!=23,3:5],2,mean)
	mnZ<-apply(M[M[,1]==23,3:5],2,mean)
	osA<-order(mnA)
	osZ<-order(mnZ)
	i1<-c(M[M[,1]!=23,2+osA[1]],M[M[,1]==23,2+osZ[1]])
	i2<-c(M[M[,1]!=23,2+osA[2]],M[M[,1]==23,2+osZ[2]])
	azdifs<-i2-i1
	## drop less than 1.5 x 25th percentile
	qs<-quantile(M[,2],probs=c(.25,.75))
	minq<-qs[1]-1.5*(qs[2]-qs[1])
	drop<-which(M[,2]<minq)
	sigs[drop]<-NA
	difs[drop]<-NA
	mdifs[drop]<-NA
	o<-cbind(ch=M[,1],sig=sigs,delta=difs,delta2=mdifs,delta3=azdifs)
	return(o)
}

wo<-vector("list")

## 10 kb windows BTBxGNPxSIN 
# A = BTB, B = GNP, C = SIN

scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
scDepth<-vector("list",23)
scID<-vector("list",23)
for(ch in 1:23){
	inf<-paste("winBTBxGNPxSIN/winout",ch,".tsv",sep="")
	dat<-read.table(inf,header=TRUE)
	L<-dim(dat)[1]
	scAB[[ch]]<-dat$A.B
	scAC[[ch]]<-dat$A.C
	scBC[[ch]]<-dat$B.C
	scDepth[[ch]]<-dat$QuartetCnt
	scID[[ch]]<-rep(ch,L)
}

## (BTB,SIN), GNP
mat<-cbind(unlist(scID),unlist(scDepth),unlist(scAC),unlist(scAB),unlist(scBC))

apply(mat[,3:5],2,mean)
wo[[1]]<-summarizeWins(mat)
##[1] 0.22553151 0.04473787 0.03596904

## 10 kb windows
# A = BCR, B = GNP, C = SIN

scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
scDepth<-vector("list",23)
scID<-vector("list",23)
for(ch in 1:23){
	inf<-paste("winBCRxGNPxSIN/winout",ch,".tsv",sep="")
	dat<-read.table(inf,header=TRUE)
	L<-dim(dat)[1]
	scAB[[ch]]<-dat$A.B
	scAC[[ch]]<-dat$A.C
	scBC[[ch]]<-dat$B.C
	scDepth[[ch]]<-dat$QuartetCnt
	scID[[ch]]<-rep(ch,L)
}

## ## (BCR,SIN), GNP
mat<-cbind(unlist(scID),unlist(scDepth),unlist(scAC),unlist(scAB),unlist(scBC))

apply(mat[,3:5],2,mean)
#[1] 0.22565410 0.04367563 0.03604101
wo[[2]]<-summarizeWins(mat)


## 10 kb windows 
## A = HNV, B = GNP, C = SIN
scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
scDepth<-vector("list",23)
scID<-vector("list",23)
for(ch in 1:23){
	inf<-paste("winHNVxGNPxSIN/winout",ch,".tsv",sep="")
	dat<-read.table(inf,header=TRUE)
	L<-dim(dat)[1]
	scAB[[ch]]<-dat$A.B
	scAC[[ch]]<-dat$A.C
	scBC[[ch]]<-dat$B.C
	scDepth[[ch]]<-dat$QuartetCnt
	scID[[ch]]<-rep(ch,L)
}

## rearrange (HNV,SIN),GNP
mat<-cbind(unlist(scID),unlist(scDepth),unlist(scAC),unlist(scAB),unlist(scBC))

apply(mat[,3:5],2,mean)
#[1] 0.10285230 0.06005708 0.05010287
wo[[3]]<-summarizeWins(mat)


## 10 kb windows 
# A = ABM, B = SIN, C = TBY
scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
scDepth<-vector("list",23)
scID<-vector("list",23)
for(ch in 1:23){
	inf<-paste("winABMxSINxTBY/winout",ch,".tsv",sep="")
	dat<-read.table(inf,header=TRUE)
	L<-dim(dat)[1]
	scAB[[ch]]<-dat$A.B
	scAC[[ch]]<-dat$A.C
	scBC[[ch]]<-dat$B.C
	scDepth[[ch]]<-dat$QuartetCnt
	scID[[ch]]<-rep(ch,L)
}

mat<-cbind(unlist(scID),unlist(scDepth),unlist(scAB),unlist(scAC),unlist(scBC))

apply(mat[,3:5],2,mean)
#[1] 0.47893530 0.05073642 0.05181248
wo[[4]]<-summarizeWins(mat)

## 10 kb windows 
## A = BHP, B = SIN, C = YG

scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
scDepth<-vector("list",23)
scID<-vector("list",23)
for(ch in 1:23){
	inf<-paste("winBHPxSINxYG/winout",ch,".tsv",sep="")
	dat<-read.table(inf,header=TRUE)
	L<-dim(dat)[1]
	scAB[[ch]]<-dat$A.B
	scAC[[ch]]<-dat$A.C
	scBC[[ch]]<-dat$B.C
	scDepth[[ch]]<-dat$QuartetCnt
	scID[[ch]]<-rep(ch,L)
}

mat<-cbind(unlist(scID),unlist(scDepth),unlist(scAB),unlist(scAC),unlist(scBC))

apply(mat[,3:5],2,mean)
#[1] 0.53677299 0.06128499 0.03884973
wo[[5]]<-summarizeWins(mat)

## 10 kb windows 
# A = BCR, B = BTB, C = YG
scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
scDepth<-vector("list",23)
scID<-vector("list",23)
for(ch in 1:23){
	inf<-paste("winBCRxBTBxYG/winout",ch,".tsv",sep="")
	dat<-read.table(inf,header=TRUE)
	L<-dim(dat)[1]
	scAB[[ch]]<-dat$A.B
	scAC[[ch]]<-dat$A.C
	scBC[[ch]]<-dat$B.C
	scDepth[[ch]]<-dat$QuartetCnt
	scID[[ch]]<-rep(ch,L)
}

mat<-cbind(unlist(scID),unlist(scDepth),unlist(scAB),unlist(scAC),unlist(scBC))

apply(mat[,3:5],2,mean)
#[1] 0.31465519 0.02925973 0.02891468
wo[[6]]<-summarizeWins(mat)

## 10 kb windows 
## A = LS, B = YG, C = SIN
scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
scDepth<-vector("list",23)
scID<-vector("list",23)
for(ch in 1:23){
	inf<-paste("winLSxYGxSIN/winout",ch,".tsv",sep="")
	dat<-read.table(inf,header=TRUE)
	L<-dim(dat)[1]
	scAB[[ch]]<-dat$A.B
	scAC[[ch]]<-dat$A.C
	scBC[[ch]]<-dat$B.C
	scDepth[[ch]]<-dat$QuartetCnt
	scID[[ch]]<-rep(ch,L)
}

mat<-cbind(unlist(scID),unlist(scDepth),unlist(scAB),unlist(scAC),unlist(scBC))

apply(mat[,3:5],2,mean)
#[1] 1.15329780 0.02080008 0.01757438
wo[[7]]<-summarizeWins(mat)

######## significance testing ###############
sigTest<-function(W=NA,thin=5,reps=1000){
	nw<-dim(W)[1]
	sigMn<-tapply(INDEX=W[,1]==23,X=W[,2],mean,na.rm=TRUE)
	difMn<-tapply(INDEX=W[,1]==23,X=W[,3],mean,na.rm=TRUE)
	dif2Mn<-tapply(INDEX=W[,1]==23,X=W[,4],mean,na.rm=TRUE)
	nullSig<-rep(NA,reps)
	nullDif<-rep(NA,reps)
	nullDif2<-rep(NA,reps)
	for(j in 1:reps){
		k<-sample(1:thin,1)
		xx<-seq(from=k,to=nw,by=thin)
		sub<-W[xx,]
		sigObs<-tapply(INDEX=sub[,1]==23,X=sub[,2],mean,na.rm=TRUE)
		difObs<-tapply(INDEX=sub[,1]==23,X=sub[,3],mean,na.rm=TRUE)
		dif2Obs<-tapply(INDEX=sub[,1]==23,X=sub[,4],mean,na.rm=TRUE)
		sub[,1]<-sample(sub[,1],dim(sub)[1],replace=FALSE)
		sigRep<-tapply(INDEX=sub[,1]==23,X=sub[,2],mean,na.rm=TRUE)
		difRep<-tapply(INDEX=sub[,1]==23,X=sub[,3],mean,na.rm=TRUE)
		dif2Rep<-tapply(INDEX=sub[,1]==23,X=sub[,4],mean,na.rm=TRUE)
		nullSig[j]<-(sigObs[2]-sigObs[1])-(sigRep[2]-sigRep[1]) 
		nullDif[j]<-(difObs[2]-difObs[1])-(difRep[2]-difRep[1]) 
		nullDif2[j]<-(dif2Obs[2]-dif2Obs[1])-(dif2Rep[2]-dif2Rep[1]) 
	}
	out<-list(sigMn=sigMn,difMn=difMn,dif2Mn=dif2Mn,nullSig=nullSig,nullDif=nullDif,
		  nullDif2=nullDif2)
	cat(mean(nullSig > 0),"\n")
	cat(mean(nullDif > 0),"\n")
	cat(mean(nullDif2 > 0),"\n")
	return(out)
}

ilsTest<-function(W=NA,thin=3,reps=1000){
	nw<-dim(W)[1]
	dif2Mn<-mean(W[,4],na.rm=TRUE)
	dif3Mn<-tapply(INDEX=W[,1]==23,X=W[,5],mean,na.rm=TRUE)
	## overall, A, Z
	difs<-c(dif2Mn,dif3Mn)
	null<-matrix(NA,nrow=reps,ncol=3)
	for(j in 1:reps){
		k<-sample(1:thin,1)
		xx<-seq(from=k,to=nw,by=thin)
		xx<-sample(xx,length(xx),replace=TRUE)
		sub<-W[xx,]
		sdif2Mn<-mean(sub[,4],na.rm=TRUE)
		sdif3Mn<-tapply(INDEX=sub[,1]==23,X=sub[,5],mean,na.rm=TRUE)
		null[j,]<-c(sdif2Mn,sdif3Mn)
	}
	cat(quantile(null[,1],probs=c(.5,.025,.975)),"\n")
	cat(quantile(null[,2],probs=c(.5,.025,.975)),"\n")
	cat(quantile(null[,3],probs=c(.5,.025,.975)),"\n")
	return(null)
}

pseq<-1:7

so<-vector("list",7)
io<-vector("list",7)
for(i in 1:7){
	so[[i]]<-sigTest(wo[[i]],thin=5)
	io[[i]]<-ilsTest(wo[[i]],thin=5)
}

## figures and tables
trs<-c("((BTB,SIN),GNP)","((BCR,SIN),GNP)","((HNV,SIN),GNP)","((ABM,SIN),TBY)","((BHP,SIN),YG)",
       "((BCR,BTB),YG)","((LS,YG),SIN)")

library(xtable)
stab<-matrix(NA,nrow=7,ncol=6)
rownames(stab)<-trs[pseq]
k<-1
for(j in pseq){
	stab[k,1:2]<-round(so[[j]][[1]],2)
	stab[k,3]<-max(mean(so[[j]][[4]] <= 0),1/1000)
	stab[k,4:5]<-round(so[[j]][[3]],2)
        p2<-mean(so[[j]][[5]] <= 0)
        p2<-max(min(p2,1-p2)*2,2/1000)# two tail
	stab[k,6]<-p2
	k<-k+1
}

xtable(stab)
#% latex table generated in R 4.5.2 by xtable 1.8-4 package
#% Wed Mar 11 14:43:40 2026
#\begin{table}[ht]
#\centering
#\begin{tabular}{rrrrrrr}
#  \hline
# & 1 & 2 & 3 & 4 & 5 & 6 \\ 
#  \hline
#((BTB,SIN),GNP) & 0.30 & 0.35 & 0.03 & 0.01 & 0.04 & 0.03 \\ 
#  ((BCR,SIN),GNP) & 0.30 & 0.35 & 0.03 & 0.01 & 0.04 & 0.08 \\ 
#  ((HNV,SIN),GNP) & 0.21 & 0.25 & 0.04 & 0.00 & 0.09 & 0.00 \\ 
#  ((ABM,SIN),TBY) & 0.55 & 0.86 & 0.00 & 0.00 & -0.01 & 0.59 \\ 
#  ((BHP,SIN),YG) & 0.62 & 1.00 & 0.00 & 0.02 & 0.01 & 0.38 \\ 
#  ((BCR,BTB),YG) & 0.38 & 0.33 & 0.92 & 0.00 & -0.00 & 0.92 \\ 
#  ((LS,YG),SIN) & 1.17 & 1.57 & 0.00 & 0.00 & 0.00 & 0.69 \\ 
#   \hline
#\end{tabular}
#\end{table}


library(scales)
cs<-alpha(c("darkgray","orange"),.5)
pdf("fig_sigWinJH.pdf",width=9,height=9)
par(mfrow=c(3,3))
par(mar=c(4,5,2.5,1.5))
cl<-1.4;ca<-1.1;ct<-1.05;cm<-1.3
k<-1
for(j in pseq){
	chi<-1+(wo[[j]][,1]==23)
	plot(wo[[j]][,3],wo[[j]][,2],pch=19,col=cs[chi],xlab="Difference in weights",ylab="Sum of weights",cex.lab=cl)
	abline(v=0,lty=2)
	obs<-round(so[[j]][[1]],2)
	p<-max(mean(so[[j]][[4]] <= 0),1/1000)
	p2<-max(mean(so[[j]][[5]] <= 0),1/1000)
	p2<-min(p2,1-p2)*2# two tail
	mtext(paste("Sum Z > A, P = ",p,sep=""),1,line=-3,adj=.9,cex=.75)
	mtext(paste("Dif != 0, P = ",p2,sep=""),1,line=-1.5,adj=.9,cex=.75)
	title(main=trs[j],cex.main=cm)
	#title(main=paste("(",LETTERS[k],") ",trs[j],sep=""),cex.main=cm)
	k<-k+1
}
plot(c(0,1),c(0,1),type='n',xlab="",ylab="",axes=FALSE)
legend(x=.15,y=.8,c("Autosomes","Z chromosome"),pch=19,col=cs,bty='n',cex=1.4)
dev.off()

pdf("fig_sigWin2JH.pdf",width=9,height=9)
par(mfrow=c(3,3))
par(mar=c(4,5,2.5,1.5))
cl<-1.4;ca<-1.1;ct<-1.05;cm<-1.3
k<-1
for(j in pseq){
	chi<-1+(wo[[j]][,1]==23)
	plot(wo[[j]][,4],wo[[j]][,2],pch=19,col=cs[chi],xlab="Difference in weights",ylab="Sum of weights",cex.lab=cl)
	abline(v=0,lty=2)
	obs<-round(so[[j]][[1]],2)
	p<-max(mean(so[[j]][[4]] <= 0),1/1000)
	p2<-max(mean(so[[j]][[6]] <= 0),1/1000)
	p2<-min(p2,1-p2)*2# two tail
	mtext(paste("Sum Z > A, P = ",p,sep=""),1,line=-3,adj=.9,cex=.75)
	mtext(paste("Dif != 0, P = ",p2,sep=""),1,line=-1.5,adj=.9,cex=.75)
	title(main=trs[j],cex.main=cm)
	#title(main=paste("(",LETTERS[k],") ",trs[j],sep=""),cex.main=cm)
	k<-k+1
}
plot(c(0,1),c(0,1),type='n',xlab="",ylab="",axes=FALSE)
legend(x=.15,y=.8,c("Autosomes","Z chromosome"),pch=19,col=cs,bty='n',cex=1.4)
dev.off()

io50<-matrix(NA,nrow=7,ncol=3)
rownames(io50)<-trs[pseq]
io05<-matrix(NA,nrow=7,ncol=3)
io95<-matrix(NA,nrow=7,ncol=3)
i<-1
for(j in pseq){
	oo<-apply(io[[j]],2,quantile,probs=c(.5,.025,.975))
	io50[i,]<-oo[1,]
	io05[i,]<-oo[2,]
	io95[i,]<-oo[3,]
	i<-i+1
}

lb<-1.05*min(io05)
ub<-1.05*max(io95)
yy<-rep(seq(1,34,5),each=3)+ rep(1:3, 7) 
yym<-matrix(rev(yy-1),nrow=3,byrow=F)
yym<-rbind(yym[3,],yym[2,],yym[1,])

pdf("fig_dotWinJH.pdf",width=7,height=7)
dotchart(t(io50),xlim=c(lb,ub),col=c("black",cs),pch=19,labels=NA,xlab="Difference in weights",cex.lab=cl,cex.axis=ca)
abline(v=0,lty=2)
segments(t(io05),yym,t(io95),yym,col=c("black",cs),lwd=1.5)
dev.off()

save(list=ls(),file="winsJH.rdat")
