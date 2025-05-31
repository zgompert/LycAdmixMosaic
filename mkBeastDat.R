## generate input file for caster 

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

SSeq<-vector("list",27)
for(i in 1:N){
	SSeq[[i]]<-vector("list",length(a1f))
}

for(i in 1:N){
	cat(i,"\n")
	out<-paste("max_chrom",chrom[i],".fasta",sep="")
	a1<-as.matrix(fread(a1f[i],header=F))
	a2<-as.matrix(fread(a2f[i],header=F))
	n<-a1+a2
	p<-a2/(a1+a2) ## non-ref
	p[n < 5]<-NA
	J<-dim(p)[2]
	L<-dim(p)[1]
	snps<-as.data.frame(fread(asnp[i],header=FALSE))
	for(j in 1:J){
		nx<-as.numeric(p[,j] > .5) + 1
		ss<-rep("N",L)
		jx<-which(is.na(p[,j])==FALSE)
		for(l in jx){
			ss[l]<-snps[l,nx[l]]
		}
		SS1<-paste(ss,collapse="")
		#SSeq[[j]][[i]]<-paste(ss,collapse="")
		cat(">",ids[j,1],"\n",file=out,append=TRUE,sep="")
		cat(SS1,"\n",file=out,append=TRUE,sep="")
	}
	
}



