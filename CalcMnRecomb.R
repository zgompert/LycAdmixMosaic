## calculate mean recombination rate per chromosome for analyses
rf<-list.files(pattern="winavg")
Rdat<-vector("list",length(rf))
for(i in 1:length(rf)){
	Rdat[[i]]<-read.table(rf[i],sep=",",header=TRUE)
}

rmat<-matrix(NA,nrow=22,ncol=length(rf))
for(i in 1:length(rf)){
	rmat[,i]<-tapply(X=Rdat[[i]][,4],INDEX=Rdat[[i]][,1],median)
}

cor(rmat)
#          [,1]      [,2]      [,3]      [,4]      [,5]
#[1,] 1.0000000 0.6026119 0.6018043 0.4718252 0.7201437
#[2,] 0.6026119 1.0000000 0.8378906 0.5448609 0.7962067
#[3,] 0.6018043 0.8378906 1.0000000 0.5970850 0.8313422
#[4,] 0.4718252 0.5448609 0.5970850 1.0000000 0.6006778
#[5,] 0.7201437 0.7962067 0.8313422 0.6006778 1.0000000

mnRho<-apply(rmat,1,mean)

cor(mnRho,rmat[,1])
#[1] 0.6881044
cor(mnRho,rmat[,2])
#[1] 0.9220997
cor(mnRho,rmat[,3])
#[1] 0.9573857
cor(mnRho,rmat[,4])
#[1] 0.6337681
cor(mnRho,rmat[,5])
#[1] 0.9339974

ssz<-read.table("../TreeMx/ChromNameSize.txt",header=FALSE)

cor.test(ssz[1:22,3],mnRho)
#	Pearson's product-moment correlation
#
#data:  ssz[1:22, 3] and mnRho
#t = -5.1772, df = 20, p-value = 4.579e-05
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.8933417 -0.4921858
#sample estimates:
#       cor 
#-0.7567594 

nms<-c("L. anna","L. idas","L. melissa","Sierra Nevada","Warner Mts.")

## median rhos by chromsome size
pdf("sfig_RhoCrhom.pdf",width=6,height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:5){
        mns<-rmat[,i]
        lns<-ssz[1:22,3]
        corr<-round(cor(lns,mns),3)
        cl<-1.5;ca<-1.05;cm<-1.4
        plot(lns,mns,pch=19,xlab="Chromosome size (bps)",ylab="Median recombination",cex.lab=cl,cex.axis=ca)
        mtext(text=paste("R = ",corr,sep=""),side=3,line=-2,adj=.82)
        title(main=paste("(",LETTERS[i],") ",nms[i],sep=""),cex.main=cm)
}
corr<-round(cor(lns,mnRho),3)
plot(lns,mnRho,pch=19,xlab="Chromosome size (bps)",ylab="Median recombination",cex.lab=cl,cex.axis=ca)
mtext(text=paste("R = ",corr,sep=""),side=3,line=-2,adj=.82)
title(main=paste("(",LETTERS[6],") Average",sep=""),cex.main=cm)
dev.off()

colnames(rmat)<-c("rhoLa","rhoLi","rhoLm","rhoSn","rhoWm")
o<-cbind(ch=1:22,mnRho,rmat)
write.table(file="ChromRho.txt",o,row.names=FALSE,col.names=TRUE,quote=FALSE)

