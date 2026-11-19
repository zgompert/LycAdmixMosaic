## summarize genetic maps

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



## pos by map
pdf("sfig_chromMaps.pdf",width=8,height=12)
par(mfrow=c(6,4))
par(mar=c(4,4,1,1))
for(i in 1:23){
	plot(as.numeric(posL[[i]])/1e6,mapL[[i]][,2],pch=19,xlab="Position (Mb)",ylab="Map position (cM)",cex.lab=1.2,cex.axis=1)
	mtext(paste("Chrom. ",i,sep=""),side=3,adj=.1,line=-2,cex=1.1)
}
dev.off()

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
#  1.615   1.844   2.154   2.470   2.891   4.163

summary(msz)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 34.83   43.90   49.38   50.19   55.13   73.76 

summary(psz)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  10.17   18.14   21.07   21.51   25.26   33.22   

cor.test(psz,msz)
#t = -0.037939, df = 21, p-value = 0.9701
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.4190509  0.4053065
#sample estimates:
#         cor 
#-0.008278777 

cor.test(psz[-23],msz[-23])
#t = -0.01544, df = 20, p-value = 0.9878
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.4244431  0.4187655
#sample estimates:
#         cor 
#-0.003452467

summary(lm(msz~psz))
#Residuals:
#     Min       1Q   Median       3Q      Max 
#-15.5249  -6.2514  -0.7617   5.0152  23.5138 

#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 50.50202    8.53498   5.917 7.13e-06 ***
#psz         -0.01467    0.38662  -0.038     0.97    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 9.238 on 21 degrees of freedom
#Multiple R-squared:  6.854e-05,	Adjusted R-squared:  -0.04755 
#F-statistic: 0.001439 on 1 and 21 DF,  p-value: 0.9701


## negative correlation between recomb rate and size
cor.test(psz,rrate)
#data:  psz and rrate
#t = -5.8054, df = 21, p-value = 9.203e-06
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.9044829 -0.5509942
#sample estimates:
#       cor 
#-0.7849237 

cor.test(psz[-23],rrate[-23])
#data:  psz[-23] and rrate[-23]
#t = -5.5053, df = 20, p-value = 2.176e-05
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.9024612 -0.5270369
#sample estimates:
#       cor 
#-0.7761765 


## compare with Rho estimates... these include Ne of course
rho<-read.table("../Recomb_iSMC/ChromRho.txt",header=TRUE)

cor.test(rrate[-23],rho$mnRho)
#t = 2.3834, df = 20, p-value = 0.02718
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.06076792 0.74433754
#sample estimates:
#      cor 
#0.4703265 

cor.test(rrate[-23],rho$rhoLm)
#data:  rrate[-23] and rho$rhoLm
#t = 2.3005, df = 20, p-value = 0.03231
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.04439307 0.73692464
#sample estimates:
#      cor 
#0.4574399 

cor.test(rho$mnRho,psz[-23])
#data:  rho$mnRho and psz[-23]
#t = -5.1765, df = 20, p-value = 4.587e-05
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.8933204 -0.4921060
#sample estimates:
#       cor 
#-0.7567143 

cor.test(rho$rhoLm,psz[-23])
#t = -4.451, df = 20, p-value = 0.0002452
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.8686795 -0.4039666
#sample estimates:
#       cor 
#-0.7054293 

pdf("sfig_rhoVr.pdf",width=4,height=8)
par(mfrow=c(2,1))
par(mar=c(4.5,4.5,1,1))
plot(rrate[-23],rho$rhoLm,type='n',xlab="Recombination (cM/Mb)",ylab="Rho (L. melissa)",cex.lab=1.2)
text(rrate[-23],rho$rhoLm,1:22)
plot(rrate[-23],rho$mnRho,type='n',xlab="Recombination (cM/Mb)",ylab="Rho (mean)",cex.lab=1.2)
text(rrate[-23],rho$mnRho,1:22)
dev.off()

