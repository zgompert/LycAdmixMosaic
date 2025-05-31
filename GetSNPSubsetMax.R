## identify SNPs from the alignments that are variable in the alignments
## and subsample these

library(data.table)

ids<-read.table("../GenData/IDs.txt",header=FALSE)
miss<-vector("list",23)
reps<-grep(pattern="rep",ids[,1])
for(i in 1:23){
	ifile<-paste("text_max_chrom",i,".fasta",sep="")
	dat<-fread(ifile,header=FALSE)
	miss[[i]]<-apply(dat[-reps,]==5,2,mean)

}
## retain 0.025 percent with the appropriate conditions
keepSNPs<-vector("list",23)
prop<-0.00025
for(i in 1:23){
	xx<-which(miss[[i]]==0) ## no missing data
	keepSNPs[[i]]<-sort(sample(xx,floor(length(miss[[i]])*prop),replace=FALSE))
}
for(i in 1:23){
	out<-paste("keepSNPs_max_chrom",i,sep="")
	write.table(keepSNPs[[i]],file=out,row.names=FALSE,col.names=FALSE,quote=FALSE)
}
save(list=ls(),file="snps_max.rdat")


