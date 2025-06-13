## identify SNPs from the alignments that are variable in the alignments
## and subsample these

library(data.table)

miss<-vector("list",23)
for(i in 1:23){
	ifile<-paste("text_max_chrom",i,".fasta",sep="")
	dat<-fread(ifile,header=FALSE)
	miss[[i]]<-apply(dat==5,2,mean)

}
## retain 0.035 percent with the appropriate conditions
keepSNPs<-vector("list",23)
prop<-0.00035
for(i in 1:23){
	xx<-which(miss[[i]]==0) ## no missing data
	keepSNPs[[i]]<-sort(sample(xx,floor(length(miss[[i]])*prop),replace=FALSE))
}
for(i in 1:23){
	out<-paste("keepSNPs_max_chrom",i,sep="")
	write.table(keepSNPs[[i]],file=out,row.names=FALSE,col.names=FALSE,quote=FALSE)
}
save(list=ls(),file="snps_max.rdat")


