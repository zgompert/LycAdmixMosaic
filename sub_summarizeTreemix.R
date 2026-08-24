source("plotting_funcs.R")
pve<-matrix(NA,nrow=23,ncol=9)
mig<-seq(0,8)
pdf("sub_covResidDump.pdf",width=6,height=6)
for(j in 1:23){
	for(i in mig){
		lf<-paste("sub_tro_ch",j,"_m",i,".llik",sep="")
		if(file.exists(lf) & file.size(lf) > 0){
		    prefix<-paste("sub_tro_ch",j,"_m",i,sep="")
		    R<-plot_resid(prefix,"poplist.txt")
		    W<-plot_cov(prefix,"poplist.txt")
		    n<-23
		    Rbar<-sum(R[upper.tri(R)])/(n*(n-1)/2)
		    Wbar<-sum(W[upper.tri(W)])/(n*(n-1)/2)
		    num<-sum((R[upper.tri(R)]-Rbar)^2)
		    dnom<-sum((W[upper.tri(W)]-Wbar)^2)
		    pve[j,i+1]<-1-(num/dnom)
		    plot_tree(prefix)
		}
	}
}
dev.off()


ll<-matrix(NA,nrow=23,ncol=9)
for(j in 1:23){
	for(i in mig){
		lf<-paste("sub_tro_ch",j,"_m",i,".llik",sep="")
		if(file.exists(lf) & file.size(lf) > 0){
			ld<-read.table(lf)
			ll[j,i+1]<-ld[2,7]
		}
	}
}
