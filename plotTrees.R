library(ape)
## ape version 5.8

pdf("casterTrees.pdf",width=9,height=9)
par(mfrow=c(3,3))
par(mar=c(1,1,3,1))
for(i in 1:23){
	inf<-paste("cout",i,sep="")
	tree<-read.tree(inf)
	plot.phylo(tree,cex=.7,use.edge.length=FALSE)
	title(main=paste("Chrom.",i),cex.main=1.3)
}

dev.off()

pdf("casterTreesMax.pdf",width=9,height=9)
par(mfrow=c(3,3))
par(mar=c(1,1,3,1))
for(i in 1:23){
	inf<-paste("cout_max_",i,sep="")
	tree<-read.tree(inf)
	plot.phylo(tree,cex=.7,use.edge.length=FALSE,type="cladogram")
	title(main=paste("Chromosome",i),cex.main=1.3)
}

dev.off()

## trying to clean things up
trees<-vector("list",23)
for(i in 1:23){
	inf<-paste("cout_max_",i,sep="")
	trees[[i]]<-read.tree(inf)
}

ref_tree<-trees[[1]]
ref_order<-ref_tree$tip.label
## rotateConstr is part of ape
alntrees<-lapply(trees, function(tr) rotateConstr(tr, ref_order))
pdf("casterTreesMax.pdf",width=9,height=9)
par(mfrow=c(3,3))
par(mar=c(1,1,3,1))
for(i in 1:23){
	plot.phylo(alntrees[[i]],cex=.7,use.edge.length=FALSE,type="cladogram")
	title(main=paste("Chromosome",i),cex.main=1.3)
}

dev.off()

## formated version of figure
pdf("fig_caster.pdf",width=8,height=10)
par(mfrow=c(4,6))
par(mar=c(1,1,1,1))
for(i in 1:23){
	plot.phylo(alntrees[[i]],cex=.7,use.edge.length=FALSE,type="cladogram")
	mtext(paste("Chr.",i),cex=1.1,line=-2,side=3,adj=.1)
}

dev.off()
