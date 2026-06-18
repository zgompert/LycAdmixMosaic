## plot window-averaged recombination landscapes, and correlations in recombination rates

afiles<-list.files(pattern="ch_anna_ismc.rho")
ifiles<-list.files(pattern="ch_idas_ismc.rho")
mfiles<-list.files(pattern="ch_melissa_ismc.rho")
sfiles<-list.files(pattern="ch_sierra_ismc.rho")
wfiles<-list.files(pattern="ch_sierra_ismc.rho") ## temp, just to make the plots
#wfiles<-list.files(pattern="ch_warner_ismc.rho")

ar<-vector("list",3)
ir<-vector("list",3)
mr<-vector("list",3)
sr<-vector("list",3)
wr<-vector("list",3)


for(i in 1:3){
	ar[[i]]<-read.table(afiles[i],header=TRUE)
	ir[[i]]<-read.table(ifiles[i],header=TRUE)
	mr[[i]]<-read.table(mfiles[i],header=TRUE)
	sr[[i]]<-read.table(sfiles[i],header=TRUE)
	wr[[i]]<-read.table(wfiles[i],header=TRUE)
}

## average rho over the 3 scales: 25k, 250k 1mb
## also sort by chrom number
winRho<-function(rr=NULL){
	avrr<-rr[[3]] ## start with 25k
	for(k in 1:dim(avrr)[1]){
		mid<-(avrr[k,2]+avrr[k,3])/2
		ch<-avrr[k,1]
		r1mi<-which(rr[[1]][,1]==ch & rr[[1]][,2] <= mid & rr[[1]][,3] >= mid)
		r250ki<-which(rr[[2]][,1]==ch & rr[[2]][,2] <= mid & rr[[2]][,3] >= mid)
		avrr[k,4]<-mean(c(avrr[k,4],rr[[1]][r1mi,4],rr[[2]][r250ki,4]))
	}
	out<-avrr[order(avrr[,1]),]
	return(out)
}

winAr<-winRho(rr=ar)
winIr<-winRho(rr=ir)
winMr<-winRho(rr=mr)
winSr<-winRho(rr=sr)
winWr<-winRho(rr=wr)

wins<-list(winAr,winIr,winMr,winSr,winWr)


pdf("sfig_Rho.pdf",width=9,height=10)
par(mfrow=c(5,1))
par(mar=c(4.5,4.5,2.5,1))
cl<-1.4;ca<-1.05;cm<-1.4
ccs<-c("antiquewhite","azure")
nms<-c("L. anna","L. idas","L. melissa","Sierra Nevada","Warner Mts.")
for(i in 1:5){
	plot(wins[[i]][,4],type='n',xlab="Chromsome",ylab="Recombination",cex.lab=cl,axes=FALSE)
	xn<-1:dim(wins[[i]])[1]
	mins<-tapply(X=xn,INDEX=wins[[i]][,1],min)
	maxs<-tapply(X=xn,INDEX=wins[[i]][,1],max)
	ymax<-1.1*max(wins[[i]][,4])
	for(k in seq(1,21,2)){
		polygon(c(mins[k],maxs[k],maxs[k],mins[k]),c(0,0,ymax,ymax),col=ccs[1],border=NA)
	}
	for(k in seq(2,22,2)){
		polygon(c(mins[k],maxs[k],maxs[k],mins[k]),c(0,0,ymax,ymax),col=ccs[2],border=NA)
	}
	lines(wins[[i]][,4],lwd=.3)
	axis(2,cex.axis=ca)
	mids<-tapply(X=1:dim(wins[[i]])[1],INDEX=wins[[i]][,1],median)
	axis(1,at=mids,1:22,cex.axis=ca)
	title(main=paste("(",LETTERS[i],") ",nms[i],sep=""),cex.main=cm)
}

dev.off()


## median rhos by chromsome size
pdf("sfig_RhoCrhom.pdf",width=6,height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:5){
	mns<-tapply(X=wins[[i]][,4],INDEX=wins[[i]][,1],median)
	lns<-tapply(X=wins[[i]][,3],INDEX=wins[[i]][,1],max)
	corr<-round(cor(lns,mns),3)
	cl<-1.5;ca<-1.05;cm<-1.4
	plot(lns,mns,pch=19,xlab="Chromosome size (bps)",ylab="Median recombination",cex.lab=cl,cex.axis=ca)
	mtext(text=paste("R = ",corr,sep=""),side=3,line=-2,adj=.82)
	title(main=paste("(",LETTERS[i],") ",nms[i],sep=""),cex.main=cm)
}
dev.off()


## correlations in recombination maps
## currently 2, 4 and 5 have an earlier first window on ch1
rmat<-as.matrix(cbind(wins[[1]][,4],wins[[2]][-1,4],wins[[3]][,4],wins[[4]][-1,4],wins[[5]][-1,4]))
crmat<-cor(rmat)

library(fields)
pdf("sfig_RhoCors.pdf",width=5.2,height=4.4)
par(mar=c(3,3,1,1))
image.plot(crmat,col = hcl.colors(12, "YlOrRd", rev = TRUE),axes=FALSE,legend.lab="correlation")
axis(1,at=seq(0,1,.25),c("La","Li","Lm","SN","WM"))
axis(2,at=seq(0,1,.25),c("La","Li","Lm","SN","WM"))
box()
dev.off()

## write out windowed averages to compare with ancestry patterns
write.table(wins[[1]],file="winavg_rho_anna.csv",sep=",",row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(wins[[2]],file="winavg_rho_idas.csv",sep=",",row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(wins[[3]],file="winavg_rho_melissa.csv",sep=",",row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(wins[[4]],file="winavg_rho_sierra.csv",sep=",",row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(wins[[5]],file="winavg_rho_warner.csv",sep=",",row.names=FALSE,col.names=TRUE,quote=FALSE)



