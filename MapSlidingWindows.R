## compare window-based estimates or r with window-based rho 

## read in position files
posL<-vector("list",23)
for(i in 1:23){
	pdat<-read.table(paste("chr_",i,".pos",sep=""),header=FALSE)
	posL[[i]]<-pdat[,1]
}

## read in map files
mapL<-vector("list",23)
for(i in 1:23){
	mapL[[i]]<-read.table(paste("chr_",i,".map",sep=""),header=FALSE)
}

## get mean recombination rate by chromosome
msz<-rep(NA,23)
psz<-rep(NA,23)
for(i in 1:23){
	psz[i]<-max(as.numeric(posL[[i]]),na.rm=TRUE)/1e6 ## in Mb
	msz[i]<-max(mapL[[i]][,2])
}
rrate<-msz/psz
summary(rrate)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  1.307   1.667   1.894   2.080   2.482   3.565 

summary(msz)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  28.18   38.84   42.80   42.65   45.41   63.16 

summary(psz)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  10.17   18.14   21.07   21.51   25.26   33.22   

## read in rho estimates from R workspace
load("../Recomb_iSMC/rhoIls.rdat")

## r in 1 Mb windows
rr<-vector("list",23)
for(k in 1:23){
	cM<-as.numeric(mapL[[k]][,2])
	pos<-as.numeric(unlist(posL[[k]]))/1e6
	ub<-max(floor(pos),na.rm=TRUE)
	rr[[k]]<-rep(NA,ub)
	for(i in 1:ub){
		a<-which(pos >= (i-1) & pos < i)
		pa<-pos[a]
		cma<-cM[a]
		rr[[k]][i]<-(max(cma,na.rm=TRUE)-min(cma,na.rm=TRUE))/(max(pa,na.rm=TRUE)-min(pa,na.rm=TRUE))
	}
}

## rho in 1 Mb windows, average and melissa... no Z, just 22 autosomes
rhoMel<-vector("list",22)
rhodat<-read.csv("../Recomb_iSMC/winavg_rho_melissa.csv")

for(k in 1:22){
	ch<-which(rhodat$chrom==k)
	pos<-((rhodat$chromStart[ch] + rhodat$chromEnd[ch])/2)/1e6 ## midpoint
	rho<-rhodat$Lyc.melissa[ch]
	ub<-max(floor(pos))
	rhoMel[[k]]<-rep(NA,ub)
	for(i in 1:ub){
		a<-which(pos >= (i-1) & pos < i)
		rhoMel[[k]][i]<-mean(rho[a],na.rm=TRUE)
	}
}

## for full comparison
## because of marker locations need to drop last rho value for ch19
## and chrom Z from recomb
x<-rhoMel
x[[19]]<-x[[19]][-18]
y<-rr
y[[23]]<-NULL
cor.test(unlist(x),unlist(y))
#t = 4.6944, df = 454, p-value = 3.546e-06
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.1258171 0.3010389
#sample estimates:
#     cor 
#0.215159 
cor.test(unlist(x),unlist(y),method="spearman")
#data:  unlist(x) and unlist(y)
#S = 12899940, p-value = 7.955e-05
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.1837062 

#Warning message:
#In cor.test.default(unlist(x), unlist(y), method = "spearman") :
#  Cannot compute exact p-value with ties
 
cor.test(unlist(x),unlist(y),method="kendall")
#
#	Kendall's rank correlation tau
#
#data:  unlist(x) and unlist(y)
#z = 3.9355, p-value = 8.303e-05
#alternative hypothesis: true tau is not equal to 0
#sample estimates:
#      tau 
#0.1247014 

cors<-matrix(NA,nrow=22,ncol=4)
for(k in 1:22){
	o<-cor.test(x[[k]],y[[k]])
	cors[k,]<-c(o$estimate,o$conf.int[1:2],o$p.value)
}

## fisher combined p
chi<--2 * sum(log(cors[,4]))
#[1] 67.62671
pchisq(chi,df=2*dim(cors)[1],lower.tail=FALSE) ## df = 44, 2k
#[1] 0.01255224


library(scales)
pdf("sfig_windowRhoR.pdf",width=4.5,height=9)
cl<-1.4;ca<-1;cm=1.22
par(mfrow=c(2,1))
par(mar=c(4,5,2.5,1.5))
plot(unlist(x),unlist(y),xlab="Rho (4Ner)",ylab="Recombination (cM/Mb)",col=alpha("black",.5),cex.lab=cl,cex.axis=ca,pch=19)
mtext("r = 0.22 (95% CI = 0.13-0.30)",side=3,line=-1.5,adj=.1,cex=1.1) 
title(main="(a) Genome-wide comparision",cex.main=cm)

plot(1:22,cors[,1],pch=19,ylim=c(-1,1),xlab="Chromosome",ylab="Pearson correlation",cex.lab=cl,cex.axis=ca)
segments(1:22,cors[,2],1:22,cors[,3])
abline(h=0,lty=3)
mtext("Fisher combined p-value = 0.013",side=3,line=-1.5,adj=.1,cex=1.1) 
title(main="(b) Chromosome-specific comparison",cex.main=cm)
dev.off()

