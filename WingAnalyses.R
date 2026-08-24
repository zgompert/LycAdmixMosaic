## ch data
chd<-read.table("../../../../TreeMx/ChromNameSize.txt",header=FALSE)

snpsGNP<-read.table("snps_tab_GNP.txt",header=FALSE)

pipGNP<-matrix(NA,nrow=dim(snpsGNP)[1],ncol=17)
for(i in 1:17){
    ff<-paste("pip_o_lyc_GNP_ph",i,"_ch0.param.txt",sep="")
    pip<-read.table(ff,header=FALSE)
    pipGNP[,i]<-pip[,2]
}

qtnGNP<-matrix(NA,nrow=23,ncol=17)
for(i in 1:17){
	qtnGNP[,i]<-tapply(pipGNP[,i],INDEX=snpsGNP[,1],sum)
}	

par(mfrow=c(3,3))
for(i in 1:9){
	plot(chd[,3],qtnGNP[,i],pch=19,col=c(rep("gray",22),"red"))
}


snpsSIN<-read.table("snps_tab_SIN.txt",header=FALSE)

pipSIN<-matrix(NA,nrow=dim(snpsSIN)[1],ncol=17)
for(i in 1:17){
    ff<-paste("pip_o_lyc_SIN_ph",i,"_ch0.param.txt",sep="")
    pip<-read.table(ff,header=FALSE)
    pipSIN[,i]<-pip[,2]
}

qtnSIN<-matrix(NA,nrow=23,ncol=17)
for(i in 1:17){
	qtnSIN[,i]<-tapply(pipSIN[,i],INDEX=snpsSIN[,1],sum)
}	

par(mfrow=c(3,3))
for(i in 1:9){
	plot(chd[,3],qtnSIN[,i],pch=19,col=c(rep("gray",22),"red"))
}

snpsYBG<-read.table("snps_tab_YBG.txt",header=FALSE)

pipYBG<-matrix(NA,nrow=dim(snpsYBG)[1],ncol=17)
for(i in 1:17){
    ff<-paste("pip_o_lyc_YBG_ph",i,"_ch0.param.txt",sep="")
    pip<-read.table(ff,header=FALSE)
    pipYBG[,i]<-pip[,2]
}

qtnYBG<-matrix(NA,nrow=23,ncol=17)
for(i in 1:17){
        qtnYBG[,i]<-tapply(pipYBG[,i],INDEX=snpsYBG[,1],sum)
}

par(mfrow=c(3,3))
for(i in 1:9){
        plot(chd[,3],qtnYBG[,i],pch=19,col=c(rep("gray",22),"red"))
}

