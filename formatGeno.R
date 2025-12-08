## remove invariant SNPs and create geno files
## female Z coded as 0/2
library(data.table)

## read in SNP data
sids<-read.table("snps.txt",header=FALSE)
lgs<-read.table("lgs.txt",header=TRUE)
lgs[24,2]<-23
L<-dim(sids)[1]
chrom<-rep(NA,L)
for(i in c(1:21,23:24)){
	a<-which(sids[,1]==lgs[i,1])
	chrom[a]<-lgs[i,2]
}
## read in autosomal genotype data
adat_GNP<-as.matrix(fread("pntest_GNP_filtered_lyc_wings.txt",header=FALSE))
adat_SIN<-as.matrix(fread("pntest_SIN_filtered_lyc_wings.txt",header=FALSE))
adat_YBG<-as.matrix(fread("pntest_YBG_filtered_lyc_wings.txt",header=FALSE))

pGNP<-apply(adat_GNP,1,mean)/2
pSIN<-apply(adat_SIN,1,mean)/2
pYBG<-apply(adat_YBG,1,mean)/2

keepGNP<-which(chrom<23 & pGNP > 0.01 & pGNP < 0.99)
keepSIN<-which(chrom<23 & pSIN > 0.01 & pSIN < 0.99)
keepYBG<-which(chrom<23 & pYBG > 0.01 & pYBG < 0.99)

length(keepGNP)
#[1] 44227
length(keepSIN)
#[1] 40725
length(keepYBG)
#[1] 11254


## read in sex chromosome genotype data
sdat_GNP<-as.matrix(fread("spntest_GNP_filtered_lyc_wings.txt",header=FALSE))
sdat_SIN<-as.matrix(fread("spntest_SIN_filtered_lyc_wings.txt",header=FALSE))
sdat_YBG<-as.matrix(fread("spntest_YBG_filtered_lyc_wings.txt",header=FALSE))

spGNP<-apply(sdat_GNP,1,mean)/2
spSIN<-apply(sdat_SIN,1,mean)/2
spYBG<-apply(sdat_YBG,1,mean)/2

skeepGNP<-which(chrom==23 & spGNP > 0.01 & spGNP < 0.99)
skeepSIN<-which(chrom==23 & spSIN > 0.01 & spSIN < 0.99)
skeepYBG<-which(chrom==23 & spYBG > 0.01 & spYBG < 0.99)

length(skeepGNP)
#[1] 2673
length(skeepSIN)
#[1] 2490
length(skeepYBG)
#[1] 584


## merge the two data sets and drop rare SNPs
gGNP<-as.matrix(rbind(adat_GNP[keepGNP,],sdat_GNP[skeepGNP,]))
#[1] 46900    98
gSIN<-as.matrix(rbind(adat_SIN[keepSIN,],sdat_SIN[skeepSIN,]))
#[1] 43215    97
gYBG<-as.matrix(rbind(adat_YBG[keepYBG,],sdat_YBG[skeepYBG,]))
#[1] 11838   100

ox<-chrom[c(keepGNP,skeepGNP)]
oGNP<-data.frame(ch=chrom[c(keepGNP,skeepGNP)][order(ox)],pos=sids[c(keepGNP,skeepGNP),2][order(ox)],A=rep("A",length(ox)),
		 G=rep("G",length(ox)),gGNP[order(ox),])
ox<-chrom[c(keepSIN,skeepSIN)]
oSIN<-data.frame(ch=chrom[c(keepSIN,skeepSIN)][order(ox)],pos=sids[c(keepSIN,skeepSIN),2][order(ox)],A=rep("A",length(ox)),
                 G=rep("G",length(ox)),gSIN[order(ox),])
ox<-chrom[c(keepYBG,skeepYBG)]
oYBG<-data.frame(ch=chrom[c(keepYBG,skeepYBG)][order(ox)],pos=sids[c(keepYBG,skeepYBG),2][order(ox)],A=rep("A",length(ox)),
                 G=rep("G",length(ox)),gYBG[order(ox),])
write.table(oGNP,file="GNP.geno",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(oSIN,file="SIN.geno",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(oYBG,file="YBG.geno",quote=FALSE,row.names=FALSE,col.names=FALSE)
## allele frequencies specifically for variance explained by chromosomes
write.table(apply(oGNP[,-c(1,2)],1,mean)/2,file="gemmaP_GNP.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(apply(oSIN[,-c(1,2)],1,mean)/2,file="gemmaP_SIN.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(apply(oYBG[,-c(1,2)],1,mean)/2,file="gemmaP_YBG.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
