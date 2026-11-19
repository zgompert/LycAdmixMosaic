## get M/F from coverage on Z vs autosomes
dp<-as.matrix(read.table("depth.txt",header=FALSE))
snps<-read.table("snpPos.txt",header=FALSE)

z<-which(snps[,1]=="Scaffold_1631;HRSCAF_2129")
a<-which(snps[,1]!="Scaffold_1631;HRSCAF_2129")
mnZ<-apply(dp[z,],2,mean)
mnA<-apply(dp[a,],2,mean)
plot(mnA,mnZ,pch=19)
abline(a=0,b=1)
hist(mnZ/mnA)

## 0.9 is a safe cutoff
plot(sort(mnZ/mnA))
N<-length(mnZ)
sex<-rep(NA,N)
sex[(mnZ/mnA) > .9]<-"M"
sex[(mnZ/mnA) < .9]<-"F"
write.table(sex,file="GeneticSex.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
