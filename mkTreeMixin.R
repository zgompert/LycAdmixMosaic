## make input for treemix
library(data.table)

a1f<-list.files(pattern="ad1_fff")
a2f<-a1f
a2f<-gsub("ad1","ad2",a2f)
asnp<-a1f
asnp<-gsub("ad1","snps",asnp)
N<-length(a1f)
ids<-read.table("IDs.txt",header=FALSE)
temp<-gsub("ad1_fff_lycSpecPool_chrom","",a1f)
chrom<-gsub(".txt","",temp)
idx<-sub(pattern="Lyc-",x=ids[,1],replacement="")
reps<-grep(pattern="rep",x=idx)

for(i in 1:N){
        cat(i,"\n")
        out<-paste("chrom",chrom[i],".fasta",sep="")
        a1<-as.matrix(fread(a1f[i],header=F))
        a2<-as.matrix(fread(a2f[i],header=F))
	L<-dim(a1)[1]
	J<-dim(a1)[2]
	## combine
	combPa<-paste(a1,a2,sep=",")
	PaMat<-matrix(combPa,nrow=L,ncol=J,byrow=FALSE)
	colnames(PaMat)<-c(idx)
	write.table(file=paste("treemix_in_ch",chrom[i],".txt",sep=""),PaMat[,-reps],quote=FALSE,row.names=FALSE)
}
