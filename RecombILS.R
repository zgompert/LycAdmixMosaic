## these analyses test for a relationship between phylogenetic signal and population recombination rates across autosomes
## this is done in general, and specifically for interchromosomal variation
## and this is done with each of the recombination maps (some of course are more relevant than others for any given set)

## read in the recombination maps
rfiles<-list.files(pattern="winavg")

Nrho<-length(rfiles)
rho<-vector("list",Nrho)
for(i in 1:Nrho){
	rho[[i]]<-read.table(rfiles[i],header=TRUE,sep=",")
}

## ancestry data for ILS, from ILS script
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
	o<-cbind(ch=M[,1],pos=M[,6],sig=sigs,delta=difs,delta2=mdifs,delta3=azdifs)
	return(o)
}

## add rho data to ancestry data
addRho<-function(M=NULL,R=rho){
	rmat<-matrix(NA,nrow=dim(M)[1],ncol=length(R))
	for(j in 1:dim(M)[1]){
		for(k in 1:length(R)){
			## 5k windows
			xx<-which(R[[k]][,1] == M[j,1] & R[[k]][,2] <= M[j,2] & R[[k]][,3] >= M[j,2])
			if(length(xx)==1){
				rmat[j,k]<-R[[k]][xx,4]
			}
		}
	}
	colnames(rmat)<-c("rA","rI","rM","rS","rW")
	out<-as.matrix(cbind(M,rmat))
	return(out)
}

wo<-vector("list")

## 10 kb windows CLHxTICxBHP
# A = CLH, B = TIC, C = BHP

scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
scDepth<-vector("list",23)
scID<-vector("list",23)
scPos<-vector("list",23)
for(ch in 1:23){
	inf<-paste("../Caster/winCLHxTICxBHP/winout",ch,".tsv",sep="")
	dat<-read.table(inf,header=TRUE)
	L<-dim(dat)[1]
	scAB[[ch]]<-dat$A.B
	scAC[[ch]]<-dat$A.C
	scBC[[ch]]<-dat$B.C
	scDepth[[ch]]<-dat$QuartetCnt
	scPos[[ch]]<-dat$pos
	scID[[ch]]<-rep(ch,L)
}

mat<-cbind(unlist(scID),unlist(scDepth),unlist(scAB),unlist(scAC),unlist(scBC),unlist(scPos))

apply(mat[,3:5],2,mean)
wo[[1]]<-summarizeWins(mat)
## drop 23
wo[[1]]<-wo[[1]][wo[[1]][,1] !=23,]

wo[[1]]<-addRho(M=wo[[1]],R=rho)


## 10 kb windows
# A = CP, B = YG, C = TIC

scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
scDepth<-vector("list",23)
scID<-vector("list",23)
scPos<-vector("list",23)
for(ch in 1:23){
	inf<-paste("../Caster/winCPxYGxTIC/winout",ch,".tsv",sep="")
	dat<-read.table(inf,header=TRUE)
	L<-dim(dat)[1]
	scAB[[ch]]<-dat$A.B
	scAC[[ch]]<-dat$A.C
	scBC[[ch]]<-dat$B.C
	scDepth[[ch]]<-dat$QuartetCnt
	scPos[[ch]]<-dat$pos
	scID[[ch]]<-rep(ch,L)
}

mat<-cbind(unlist(scID),unlist(scDepth),unlist(scAC),unlist(scAB),unlist(scBC),unlist(scPos))

apply(mat[,3:5],2,mean)
wo[[2]]<-summarizeWins(mat)
## drop 23
wo[[2]]<-wo[[2]][wo[[2]][,1] !=23,]

wo[[2]]<-addRho(M=wo[[2]],R=rho)


## 10 kb windows
# A = MR, B = CP, C = YG

scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
scDepth<-vector("list",23)
scID<-vector("list",23)
scPos<-vector("list",23)
for(ch in 1:23){
	inf<-paste("../Caster/winMRxCPxYG/winout",ch,".tsv",sep="")
	dat<-read.table(inf,header=TRUE)
	L<-dim(dat)[1]
	scAB[[ch]]<-dat$A.B
	scAC[[ch]]<-dat$A.C
	scBC[[ch]]<-dat$B.C
	scDepth[[ch]]<-dat$QuartetCnt
	scPos[[ch]]<-dat$pos
	scID[[ch]]<-rep(ch,L)
}

mat<-cbind(unlist(scID),unlist(scDepth),unlist(scAC),unlist(scAB),unlist(scBC),unlist(scPos))

apply(mat[,3:5],2,mean)
wo[[3]]<-summarizeWins(mat)
## drop 23
wo[[3]]<-wo[[3]][wo[[3]][,1] !=23,]

wo[[3]]<-addRho(M=wo[[3]],R=rho)

## 10 kb windows
# A = CLH, B = YG, C = BHP

scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
scDepth<-vector("list",23)
scID<-vector("list",23)
scPos<-vector("list",23)
for(ch in 1:23){
	inf<-paste("../Caster/winCLHxYGxBHP/winout",ch,".tsv",sep="")
	dat<-read.table(inf,header=TRUE)
	L<-dim(dat)[1]
	scAB[[ch]]<-dat$A.B
	scAC[[ch]]<-dat$A.C
	scBC[[ch]]<-dat$B.C
	scDepth[[ch]]<-dat$QuartetCnt
	scPos[[ch]]<-dat$pos
	scID[[ch]]<-rep(ch,L)
}

mat<-cbind(unlist(scID),unlist(scDepth),unlist(scAB),unlist(scAC),unlist(scBC),unlist(scPos))

apply(mat[,3:5],2,mean)
wo[[4]]<-summarizeWins(mat)
## drop 23
wo[[4]]<-wo[[4]][wo[[4]][,1] !=23,]

wo[[4]]<-addRho(M=wo[[4]],R=rho)

## 10 kb windows
# A = TIC, B = YG, C = BHP

scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
scDepth<-vector("list",23)
scID<-vector("list",23)
scPos<-vector("list",23)
for(ch in 1:23){
	inf<-paste("../Caster/winTICxYGxBHP/winout",ch,".tsv",sep="")
	dat<-read.table(inf,header=TRUE)
	L<-dim(dat)[1]
	scAB[[ch]]<-dat$A.B
	scAC[[ch]]<-dat$A.C
	scBC[[ch]]<-dat$B.C
	scDepth[[ch]]<-dat$QuartetCnt
	scPos[[ch]]<-dat$pos
	scID[[ch]]<-rep(ch,L)
}

mat<-cbind(unlist(scID),unlist(scDepth),unlist(scAB),unlist(scAC),unlist(scBC),unlist(scPos))

apply(mat[,3:5],2,mean)
wo[[5]]<-summarizeWins(mat)
## drop 23
wo[[5]]<-wo[[5]][wo[[5]][,1] !=23,]

wo[[5]]<-addRho(M=wo[[5]],R=rho)

## 10 kb windows
# A = CP, B = YG, C = BHP

scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
scDepth<-vector("list",23)
scID<-vector("list",23)
scPos<-vector("list",23)
for(ch in 1:23){
	inf<-paste("../Caster/winCPxYGxBHP/winout",ch,".tsv",sep="")
	dat<-read.table(inf,header=TRUE)
	L<-dim(dat)[1]
	scAB[[ch]]<-dat$A.B
	scAC[[ch]]<-dat$A.C
	scBC[[ch]]<-dat$B.C
	scDepth[[ch]]<-dat$QuartetCnt
	scPos[[ch]]<-dat$pos
	scID[[ch]]<-rep(ch,L)
}

mat<-cbind(unlist(scID),unlist(scDepth),unlist(scAB),unlist(scAC),unlist(scBC),unlist(scPos))

apply(mat[,3:5],2,mean)
wo[[6]]<-summarizeWins(mat)
## drop 23
wo[[6]]<-wo[[6]][wo[[6]][,1] !=23,]

wo[[6]]<-addRho(M=wo[[6]],R=rho)

## 10 kb windows
# A = EP, B = SHC, C = BHP

scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
scDepth<-vector("list",23)
scID<-vector("list",23)
scPos<-vector("list",23)
for(ch in 1:23){
	inf<-paste("../Caster/winEPxSHCxBHP/winout",ch,".tsv",sep="")
	dat<-read.table(inf,header=TRUE)
	L<-dim(dat)[1]
	scAB[[ch]]<-dat$A.B
	scAC[[ch]]<-dat$A.C
	scBC[[ch]]<-dat$B.C
	scDepth[[ch]]<-dat$QuartetCnt
	scPos[[ch]]<-dat$pos
	scID[[ch]]<-rep(ch,L)
}

mat<-cbind(unlist(scID),unlist(scDepth),unlist(scAB),unlist(scAC),unlist(scBC),unlist(scPos))

apply(mat[,3:5],2,mean)
wo[[7]]<-summarizeWins(mat)
## drop 23
wo[[7]]<-wo[[7]][wo[[7]][,1] !=23,]

wo[[7]]<-addRho(M=wo[[7]],R=rho)



