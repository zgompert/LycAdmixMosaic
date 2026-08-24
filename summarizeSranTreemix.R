source("plotting_funcs.R")
mig<-seq(0,8)
pve<-vector("list",10)
for(k in 9:10){## loop replicates
	pve[[k]]<-matrix(NA,nrow=23,ncol=9)
	pdf("SummSranDump.pdf",width=6,height=6)
	for(j in 1:23){
		for(i in mig){
		    prefix<-paste("ssub",k-1,"/sub_tro_ch",j,".gz_",i,sep="")
		    R<-plot_resid(prefix,"poplist.txt")
		    W<-plot_cov(prefix,"poplist.txt")
		    rr<-R[upper.tri(R)]
	   	    n<-23
	    	    Rbar<-sum(R[upper.tri(R)])/(n*(n-1)/2)
	            Wbar<-sum(W[upper.tri(W)])/(n*(n-1)/2)
	            num<-sum((R[upper.tri(R)]-Rbar)^2)
	            dnom<-sum((W[upper.tri(W)]-Wbar)^2)
	    	    pve[[k]][j,i+1]<-1-(num/dnom)

	        }
	}
	dev.off()

}

pve_sran<-pve

save(list=ls(),file="tmxiSran.rdat")

## load original results
load(file="tmx.rdat")


## get means for each chromosome and m
pve_sran_mns<-matrix(NA,nrow=23,ncol=9)
for(j in 1:23){for(i in 1:9){
	pve_sran_mns[j,i]<-median(c(pve_sran[[1]][j,i],pve_sran[[2]][j,i],pve_sran[[3]][j,i],
				    pve_sran[[4]][j,i],pve_sran[[5]][j,i],pve_sran[[6]][j,i],
				    pve_sran[[7]][j,i],pve_sran[[8]][j,i],pve_sran[[9]][j,i],
				    pve_sran[[10]][j,i]),na.rm=TRUE)
}}	

cors<-matrix(NA,nrow=9,ncol=2)
cors_noZ<-matrix(NA,nrow=9,ncol=2)
for(i in 1:9){
	o<-cor.test(pve[,i],pve_sran_mns[,i])
	cors[i,]<-c(o$estimate,o$p.value)
	o<-cor.test(pve[-23,i],pve_sran_mns[-23,i])
	cors_noZ[i,]<-c(o$estimate,o$p.value)
}
