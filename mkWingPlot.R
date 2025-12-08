pGNP<-read.table("../../geno/gemmaP_GNP.txt",header=FALSE)

vaGNP<-matrix(NA,nrow=17,ncol=2)
nvaGNP<-matrix(NA,nrow=17,ncol=2)
for(i in 1:17){
    ff<-paste("mav_o_lyc_GNP_ph",i,"_ch0.param.txt",sep="")
    mav<-read.table(ff,header=FALSE)
    va<-2*pGNP[,1]*(1-pGNP[,1])*mav[,2]^2
    sex<-grep(mav[,1],pattern="23:")
    vaGNP[i,2]<-sum(va[sex])
    vaGNP[i,1]<-sum(va[-sex])
    nvaGNP[i,2]<-sum(va[sex])/sum(va)
    nvaGNP[i,1]<-sum(va[-sex])/sum(va)  
}

pSIN<-read.table("../../geno/gemmaP_SIN.txt",header=FALSE)

vaSIN<-matrix(NA,nrow=17,ncol=2)
nvaSIN<-matrix(NA,nrow=17,ncol=2)
for(i in 1:17){
    ff<-paste("mav_o_lyc_SIN_ph",i,"_ch0.param.txt",sep="")
    mav<-read.table(ff,header=FALSE)
    va<-2*pSIN[,1]*(1-pSIN[,1])*mav[,2]^2
    sex<-grep(mav[,1],pattern="23:")
    vaSIN[i,2]<-sum(va[sex])
    vaSIN[i,1]<-sum(va[-sex])
    nvaSIN[i,2]<-sum(va[sex])/sum(va)
    nvaSIN[i,1]<-sum(va[-sex])/sum(va)  
}

pYBG<-read.table("../../geno/gemmaP_YBG.txt",header=FALSE)

vaYBG<-matrix(NA,nrow=17,ncol=2)
nvaYBG<-matrix(NA,nrow=17,ncol=2)
for(i in 1:17){
    ff<-paste("mav_o_lyc_YBG_ph",i,"_ch0.param.txt",sep="")
    mav<-read.table(ff,header=FALSE)
    va<-2*pYBG[,1]*(1-pYBG[,1])*mav[,2]^2
    sex<-grep(mav[,1],pattern="23:")
    vaYBG[i,2]<-sum(va[sex])
    vaYBG[i,1]<-sum(va[-sex])
    nvaYBG[i,2]<-sum(va[sex])/sum(va)
    nvaYBG[i,1]<-sum(va[-sex])/sum(va)  
}


wings<-read.table("../../../traits/resid-sizeANDcoord-15ix17-pops-NoNA.csv",sep=",",header=TRUE)
keep<-wings$pop %in% c("GNP","SIN","YBG")
wmat<-wings[keep,3:19]
o<-prcomp(wmat,center=TRUE,scale=TRUE)
#                 PC1        PC2          PC3         PC4         PC5
#a1_area    0.2295378  0.2760330 -0.095039095  0.28535496 -0.14238873
#a2_area    0.2356893  0.3500221  0.033915074  0.09810039 -0.13465150
#a3_area    0.2360654  0.3627351  0.017324161 -0.07077971  0.01955972
#a4_area    0.2455058  0.3342941 -0.001438116 -0.11124426 -0.02050613
#a5_area    0.2502273  0.3183096  0.038040477 -0.12669663  0.09927635
#a6_area    0.2345649  0.3320741 -0.095601954 -0.09009495  0.17424665
#cu23_area  0.2258990 -0.1813523 -0.313454832 -0.54148210 -0.16476008
#m_area     0.2319975 -0.2035746 -0.372596405 -0.25338092  0.26095022
#sc3_area   0.2311237 -0.1595013 -0.535657430  0.20376827  0.22318267
#sc_area    0.2401821 -0.1908011 -0.198517745  0.56616932 -0.00172881
#rs_area    0.2605326 -0.1302704  0.150464361  0.29978806 -0.06520867
#m1_area    0.2636059 -0.1616735  0.241489554  0.01746444 -0.27437394
#m2_area    0.2674027 -0.1641634  0.203885056 -0.08815428 -0.14397256
#m3_area    0.2294396 -0.1432070  0.498579057 -0.05264172  0.53069767
#cu1_area   0.2470507 -0.2051301  0.150482529  0.02400968  0.38801315
#cu21a_area 0.2548673 -0.1967992  0.145470912 -0.21342285 -0.26912971

## all individuals
wmata<-wings[,3:19]
oa<-prcomp(wmata,center=TRUE,scale=TRUE)
#                 PC1        PC2          PC3         PC4          PC5
#a1_area    0.2223730  0.3041390  0.080734751  0.06041318 -0.008935875
#a2_area    0.2352985  0.3369592 -0.014049874 -0.01287079 -0.009052393
#a3_area    0.2395017  0.3336522 -0.006174008 -0.03167964  0.038608257
#a4_area    0.2363523  0.3440337  0.001317480 -0.01644808  0.019493425
#a5_area    0.2345835  0.3447089 -0.039193743  0.01526881 -0.028064463
#a6_area    0.2229261  0.3465442  0.035688673  0.01677278 -0.017946852
#cu23_area  0.2222536 -0.1530761  0.457747632 -0.59016169 -0.300830499
#m_area     0.2457002 -0.1756415  0.310484469 -0.08594532 -0.257566300
#sc3_area   0.2295372 -0.1717520  0.448472438  0.54873523 -0.123895807
#sc_area    0.2479884 -0.2018253  0.142314578  0.45216406  0.047985851
#rs_area    0.2637466 -0.1371663 -0.263370825  0.15960368 -0.116249845
#m1_area    0.2586756 -0.1981714 -0.279923581  0.01971149  0.028285622
#m2_area    0.2660786 -0.1679997 -0.234488544  0.03550120  0.063186660
#m3_area    0.2419255 -0.1507212 -0.452235671 -0.06323807 -0.248239912
#cu1_area   0.2581789 -0.1876975 -0.159193927 -0.13295716 -0.068114712
#cu21a_area 0.2615879 -0.1878908 -0.046905279 -0.25130406  0.144879495
#X2a_area   0.2289091 -0.1458863  0.168386922 -0.14561342  0.846644149




library(scales)

pids<-as.numeric(as.factor(wings$pop[keep]))
cs<-c("#7fff00ff","#00c5cdff","#cd8500ff")

mn1<-tapply(X=o$x[,1],INDEX=pids,mean)
mn2<-tapply(X=o$x[,2],INDEX=pids,mean)

pdf("fig_wings.pdf",width=8,height=11)
layout(matrix(c(1,2,3,3,4,4,5,5),nrow=4,ncol=2,byrow=TRUE),widths=c(4,4),heights=c(4,2.33,2.33,2.33))
nms<-c("a1","a2","a3","a4","a5","a6","cu23","m","sc3","sc","rs","m1","m2","m3","cu1","cu21a","2a")
cl<-1.4;ca<-1.1
par(mar=c(4.5,4.5,2.5,.5))
plot(1:10,type='n',axes=FALSE,xlab="",ylab="")
title(main="(a) Wing pattern elements",cex.main=cl)
plot(o$x[,1],o$x[,2],xlab="PC1 (65.1%)",ylab="PC2 (13.5%)",pch=19,col=alpha(cs[pids],.7),cex.lab=cl,cex.axis=ca)
points(mn1,mn2,cex=3,pch=19,col=cs)
legend(4.5,4,c("GNP","SIN","YBG"),pch=19,col=cs,bty='n')
title(main="(b) Wing pattern PCA",cex.main=cl)
barplot(t(nvaGNP[,c(2,1)]),col=c(alpha("orange",.7),"white"),border="orange3",names.arg=nms,xlab="Pattern element",ylab="Proportion of Va",cex.lab=cl,cex.axis=ca)
title(main="(c) Pattern genetic variance in GNP",cex.main=cl)
barplot(t(nvaSIN[,c(2,1)]),col=c(alpha("orange",.7),"white"),border="orange3",names.arg=nms,xlab="Pattern element",ylab="Proportion of Va",cex.lab=cl,cex.axis=ca)
title(main="(d) Pattern genetic variance in SIN",cex.main=cl)
barplot(t(nvaYBG[,c(2,1)]),col=c(alpha("orange",.7),"white"),border="orange3",names.arg=nms,xlab="Pattern element",ylab="Proportion of Va",cex.lab=cl,cex.axis=ca)
title(main="(e) Pattern genetic variance in  YBG",cex.main=cl)
dev.off()

pidsa<-as.numeric(as.factor(wings$pop))
cs<-c("lightgray","#7fff00ff","#00c5cdff","#cd8500ff")

mn1<-tapply(X=oa$x[,1],INDEX=pidsa,mean)[-1]
mn2<-tapply(X=oa$x[,2],INDEX=pidsa,mean)[-1]

pdf("fig_wings2.pdf",width=8,height=11)
layout(matrix(c(1,2,3,3,4,4,5,5),nrow=4,ncol=2,byrow=TRUE),widths=c(4,4),heights=c(4,2.33,2.33,2.33))
nms<-c("a1","a2","a3","a4","a5","a6","cu23","m","sc3","sc","rs","m1","m2","m3","cu1","cu21a","2a")
cl<-1.4;ca<-1.1
par(mar=c(4.5,4.5,2.5,.5))
plot(1:10,type='n',axes=FALSE,xlab="",ylab="")
title(main="(a) Wing pattern elements",cex.main=cl)
plot(oa$x[,1],oa$x[,2],xlab="PC1 (65.1%)",ylab="PC2 (13.5%)",pch=19,col=alpha(cs[pidsa],.7),cex.lab=cl,cex.axis=ca)
points(mn1,mn2,cex=3,pch=19,col=cs[-1])
legend(6.5,6,c("GNP","SIN","YBG"),pch=19,col=cs[-1],bty='n')
title(main="(b) Wing pattern PCA",cex.main=cl)
barplot(t(nvaGNP[,c(2,1)]),col=c(alpha("orange",.7),"white"),border="orange3",names.arg=nms,xlab="Pattern element",ylab="Proportion of Va",cex.lab=cl,cex.axis=ca)
title(main="(c) Pattern genetic variance in GNP",cex.main=cl)
barplot(t(nvaSIN[,c(2,1)]),col=c(alpha("orange",.7),"white"),border="orange3",names.arg=nms,xlab="Pattern element",ylab="Proportion of Va",cex.lab=cl,cex.axis=ca)
title(main="(d) Pattern genetic variance in SIN",cex.main=cl)
barplot(t(nvaYBG[,c(2,1)]),col=c(alpha("orange",.7),"white"),border="orange3",names.arg=nms,xlab="Pattern element",ylab="Proportion of Va",cex.lab=cl,cex.axis=ca)
title(main="(e) Pattern genetic variance in  YBG",cex.main=cl)
dev.off()

## what about number of SNPs or chromosome sizes

mean(nvaGNP[,2]>1/23)
#[1] 0.8823529
mean(nvaSIN[,2]>1/23)
#[1] 0.8235294
mean(nvaYBG[,2]>1/23)
#[1] 0.7058824

mean(nvaGNP[,2])/(1/23)
#[1] 3.376278
mean(nvaSIN[,2])/(1/23)
#[1] 1.458684
mean(nvaYBG[,2])/(1/23)
#[1] 1.993521


## pve for si
pveGNP<-read.table("post_GNP.txt",header=FALSE)
pveSIN<-read.table("post_SIN.txt",header=FALSE)
pveYBG<-read.table("post_YBG.txt",header=FALSE)

pdf("sfit_pve.pdf",width=8,height=7)
par(mfrow=c(1,3))
par(mar=c(4.5,4,2.5,.5))
cl<-1.4;ca<-1.1
dotchart(pveGNP[,1],xlim=c(0,1),pch=19,labels=nms,cex.lab=cl,cex.axis=ca,xlab="PVE")
segments(pveGNP[,2],1:17,pveGNP[,3],1:17)
title(main="(a) Pattern PVE GNP",cex.main=cl)
dotchart(pveSIN[,1],xlim=c(0,1),pch=19,labels=nms,cex.lab=cl,cex.axis=ca,xlab="PVE")
segments(pveSIN[,2],1:17,pveSIN[,3],1:17)
title(main="(b) Pattern PVE SIN",cex.main=cl)
dotchart(pveYBG[,1],xlim=c(0,1),pch=19,labels=nms,cex.lab=cl,cex.axis=ca,xlab="PVE")
segments(pveYBG[,2],1:17,pveYBG[,3],1:17)
title(main="(c) Pattern PVE YBG",cex.main=cl)
dev.off()
