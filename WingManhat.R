## Manhattan PIP plots

minP<-0.01
nms<-c("a1","a2","a3","a4","a5","a6","cu23","m","sc3","sc","rs","m1","m2","m3","cu1","cu21a","2a")
grp1<-1:6 ## a1-a6
grp3<-7:9 ## sc3, m and cu23
grp2<-10:17 ## submarginal spots

############## GNP ###################
nGNP<-72957
pipGNP<-matrix(NA,nrow=17,ncol=nGNP)
for(i in 1:17){
	ff<-paste("pip_o_lyc_GNP_ph",i,"_ch0.param.txt",sep="")
	pip<-read.table(ff,header=FALSE)
    	pipGNP[i,]<-pip[,2]
}

chGNP<-as.numeric(gsub("(\\d+)_\\d+", "\\1", pip[,1]))

pipGNP<-pipGNP[,order(chGNP)]
chGNP<-chGNP[order(chGNP)]

xGNP<-1:nGNP
midGNP<-tapply(X=xGNP,chGNP,mean)
cs<-rep("darkgray",nGNP)
cs[chGNP==23]<-"black"

pdf("sfig_manhat_GNP.pdf",width=9,height=8)
par(mfrow=c(3,1))
par(mar=c(4.5,4.5,2.5,.5))
cl<-1.5;ca<-1.05;cm<-1.5
i<-1
mx<-max(as.numeric(pipGNP[grp1,]))

a<-which(pipGNP[i,] >= 0.01)
plot(xGNP[a],pipGNP[i,a],pch=i+1,col=cs[a],ylim=c(0,mx),xlim=c(0,nGNP),type='n',ylab="PIP",cex.lab=cl,cex.axis=ca,axes=FALSE,xlab="Chromosome")
for(i in grp1){
	a<-which(pipGNP[i,] >= 0.01)
	points(xGNP[a],pipGNP[i,a],pch=i+1,col=cs[a])
}
axis(1,at=midGNP,c(1:22,"Z"),cex.axis=ca)
axis(2,cex.axis=ca)
box()
legend(1500,.85*mx,nms[grp1],pch=grp1+1,bty='n',ncol=1)
title(main="(A)",cex.main=cm)

i<-10
mx<-max(as.numeric(pipGNP[grp2,]))

a<-which(pipGNP[i,] >= 0.01)
plot(xGNP[a],pipGNP[i,a],pch=i+1,col=cs[a],ylim=c(0,mx),xlim=c(0,nGNP),type='n',ylab="PIP",cex.lab=cl,cex.axis=ca,axes=FALSE,xlab="Chromosome")
for(i in grp2){
	a<-which(pipGNP[i,] >= 0.01)
	points(xGNP[a],pipGNP[i,a],pch=i+1,col=cs[a])
}
axis(1,at=midGNP,c(1:22,"Z"),cex.axis=ca)
axis(2,cex.axis=ca)
box()
legend(1500,.95*mx,nms[grp2],pch=grp2+1,bty='n',ncol=1)
title(main="(B)",cex.main=cm)

i<-7
mx<-max(as.numeric(pipGNP[grp3,]))

a<-which(pipGNP[i,] >= 0.01)
plot(xGNP[a],pipGNP[i,a],pch=i+1,col=cs[a],ylim=c(0,mx),xlim=c(0,nGNP),type='n',ylab="PIP",cex.lab=cl,cex.axis=ca,axes=FALSE,xlab="Chromosome")
for(i in grp3){
	a<-which(pipGNP[i,] >= 0.01)
	points(xGNP[a],pipGNP[i,a],pch=i+1,col=cs[a])
}
axis(1,at=midGNP,c(1:22,"Z"),cex.axis=ca)
axis(2,cex.axis=ca)
box()
legend(1500,.85*mx,nms[grp3],pch=grp3+1,bty='n',ncol=1)
title(main="(C)",cex.main=cm)

dev.off()

############## SIN ###################
nSIN<-66712
pipSIN<-matrix(NA,nrow=17,ncol=nSIN)
for(i in 1:17){
	ff<-paste("pip_o_lyc_SIN_ph",i,"_ch0.param.txt",sep="")
	pip<-read.table(ff,header=FALSE)
    	pipSIN[i,]<-pip[,2]
}

chSIN<-as.numeric(gsub("(\\d+)_\\d+", "\\1", pip[,1]))

pipSIN<-pipSIN[,order(chSIN)]
chSIN<-chSIN[order(chSIN)]

xSIN<-1:nSIN
midSIN<-tapply(X=xSIN,chSIN,mean)
cs<-rep("darkgray",nSIN)
cs[chSIN==23]<-"black"

pdf("sfig_manhat_SIN.pdf",width=9,height=8)
par(mfrow=c(3,1))
par(mar=c(4.5,4.5,2.5,.5))
i<-1
mx<-max(as.numeric(pipSIN[grp1,]))

a<-which(pipSIN[i,] >= 0.01)
plot(xSIN[a],pipSIN[i,a],pch=i+1,col=cs[a],ylim=c(0,mx),xlim=c(0,nSIN),type='n',ylab="PIP",cex.lab=cl,cex.axis=ca,axes=FALSE,xlab="Chromosome")
for(i in grp1){
	a<-which(pipSIN[i,] >= 0.01)
	points(xSIN[a],pipSIN[i,a],pch=i+1,col=cs[a])
}
axis(1,at=midSIN,c(1:22,"Z"),cex.axis=ca)
axis(2,cex.axis=ca)
box()
legend(1500,.85*mx,nms[grp1],pch=grp1+1,bty='n',ncol=1)
title(main="(A)",cex.main=cm)

i<-10
mx<-max(as.numeric(pipSIN[grp2,]))

a<-which(pipSIN[i,] >= 0.01)
plot(xSIN[a],pipSIN[i,a],pch=i+1,col=cs[a],ylim=c(0,mx),xlim=c(0,nSIN),type='n',ylab="PIP",cex.lab=cl,cex.axis=ca,axes=FALSE,xlab="Chromosome")
for(i in grp2){
	a<-which(pipSIN[i,] >= 0.01)
	points(xSIN[a],pipSIN[i,a],pch=i+1,col=cs[a])
}
axis(1,at=midSIN,c(1:22,"Z"),cex.axis=ca)
axis(2,cex.axis=ca)
box()
legend(1500,.95*mx,nms[grp2],pch=grp2+1,bty='n',ncol=1)
title(main="(B)",cex.main=cm)

i<-7
mx<-max(as.numeric(pipSIN[grp3,]))

a<-which(pipSIN[i,] >= 0.01)
plot(xSIN[a],pipSIN[i,a],pch=i+1,col=cs[a],ylim=c(0,mx),xlim=c(0,nSIN),type='n',ylab="PIP",cex.lab=cl,cex.axis=ca,axes=FALSE,xlab="Chromosome")
for(i in grp3){
	a<-which(pipSIN[i,] >= 0.01)
	points(xSIN[a],pipSIN[i,a],pch=i+1,col=cs[a])
}
axis(1,at=midSIN,c(1:22,"Z"),cex.axis=ca)
axis(2,cex.axis=ca)
box()
legend(1500,.85*mx,nms[grp3],pch=grp3+1,bty='n',ncol=1)
title(main="(C)",cex.main=cm)

dev.off()

############## YBG ###################
nYBG<-18161
pipYBG<-matrix(NA,nrow=17,ncol=nYBG)
for(i in 1:17){
	ff<-paste("pip_o_lyc_YBG_ph",i,"_ch0.param.txt",sep="")
	pip<-read.table(ff,header=FALSE)
    	pipYBG[i,]<-pip[,2]
}

chYBG<-as.numeric(gsub("(\\d+)_\\d+", "\\1", pip[,1]))

pipYBG<-pipYBG[,order(chYBG)]
chYBG<-chYBG[order(chYBG)]

xYBG<-1:nYBG
midYBG<-tapply(X=xYBG,chYBG,mean)
cs<-rep("darkgray",nYBG)
cs[chYBG==23]<-"black"

pdf("sfig_manhat_YBG.pdf",width=9,height=8)
par(mfrow=c(3,1))
par(mar=c(4.5,4.5,2.5,.5))
i<-1
mx<-max(as.numeric(pipYBG[grp1,]))

a<-which(pipYBG[i,] >= 0.01)
plot(xYBG[a],pipYBG[i,a],pch=i+1,col=cs[a],ylim=c(0,mx),xlim=c(0,nYBG),type='n',ylab="PIP",cex.lab=cl,cex.axis=ca,axes=FALSE,xlab="Chromosome")
for(i in grp1){
	a<-which(pipYBG[i,] >= 0.01)
	points(xYBG[a],pipYBG[i,a],pch=i+1,col=cs[a])
}
axis(1,at=midYBG,c(1:22,"Z"),cex.axis=ca)
axis(2,cex.axis=ca)
box()
legend(1500,.85*mx,nms[grp1],pch=grp1+1,bty='n',ncol=1)
title(main="(A)",cex.main=cm)

i<-10
mx<-max(as.numeric(pipYBG[grp2,]))

a<-which(pipYBG[i,] >= 0.01)
plot(xYBG[a],pipYBG[i,a],pch=i+1,col=cs[a],ylim=c(0,mx),xlim=c(0,nYBG),type='n',ylab="PIP",cex.lab=cl,cex.axis=ca,axes=FALSE,xlab="Chromosome")
for(i in grp2){
	a<-which(pipYBG[i,] >= 0.01)
	points(xYBG[a],pipYBG[i,a],pch=i+1,col=cs[a])
}
axis(1,at=midYBG,c(1:22,"Z"),cex.axis=ca)
axis(2,cex.axis=ca)
box()
legend(1500,.95*mx,nms[grp2],pch=grp2+1,bty='n',ncol=1)
title(main="(B)",cex.main=cm)

i<-7
mx<-max(as.numeric(pipYBG[grp3,]))

a<-which(pipYBG[i,] >= 0.01)
plot(xYBG[a],pipYBG[i,a],pch=i+1,col=cs[a],ylim=c(0,mx),xlim=c(0,nYBG),type='n',ylab="PIP",cex.lab=cl,cex.axis=ca,axes=FALSE,xlab="Chromosome")
for(i in grp3){
	a<-which(pipYBG[i,] >= 0.01)
	points(xYBG[a],pipYBG[i,a],pch=i+1,col=cs[a])
}
axis(1,at=midYBG,c(1:22,"Z"),cex.axis=ca)
axis(2,cex.axis=ca)
box()
legend(1500,.85*mx,nms[grp3],pch=grp3+1,bty='n',ncol=1)
title(main="(C)",cex.main=cm)

dev.off()

