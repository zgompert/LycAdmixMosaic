## identifies unique m=0 trees
library(phytools)

trf<-list.files(pattern="^tro_ch[0-9]+_m0.treeout.gz")
## list of files, in alphanumeric, 1, 10, 11, etc
N<-length(trf)
trL<-vector("list",N)
for(i in 1:N){
	trL[[i]]<-root(read.tree(trf[i]),outgroup="MEN12")
}

utr<-unique.multiPhylo(x=trL,use.edge.length=FALSE,use.tip.label=TRUE)
