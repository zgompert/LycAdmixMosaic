## read treemix infiles, randomly sample 300,000 SNPs for analysis

ff<-list.files(pattern="treemix_in")
ff<-ff[-c(14:15)] ## drop 22.1 and 22.2

for(k in 1:10){
	for(i in 1:23){
		cat("Working on",ff[i],"\n")
		dat<-read.table(ff[i],header=TRUE)
		L<-dim(dat)[1]
		snps<-sample(1:L,300000,replace=FALSE)
		do<-dat[sort(snps),]
		write.table(do,file=paste("rsub",k-1,"/sub_",ff[i],sep=""),row.names=FALSE,col.names=TRUE,quote=FALSE)
	}
}

## read treemix infiles, randomly sample start points for 
## 300,000 SNPs for analysis


for(k in 1:10){
	for(i in 1:23){
		cat("Working on",ff[i],"\n")
		dat<-read.table(ff[i],header=TRUE)
		L<-dim(dat)[1]
		st<-sample(1:(L-300000),1)
		snps<-st:(st+300000)
		do<-dat[sort(snps),]
		write.table(do,file=paste("ssub",k-1,"/sub_",ff[i],sep=""),row.names=FALSE,col.names=TRUE,quote=FALSE)
	}
}
