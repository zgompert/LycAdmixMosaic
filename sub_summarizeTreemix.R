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



## NEED TO UPDATE READ IN SNPS NUM
## pve vs chromosome size (approx, based on number of snps), excluding Z
## smaller chromosomes are better fit by bifurcating tree
cor.test(snum[-23],pve[-23,1])
#data:  snum[-23] and pve[-23, 1]
#t = -3.2901, df = 20, p-value = 0.003658
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.8114675 -0.2279422
#sample estimates:
#       cor
#-0.5926004

plot(snum[-23],pve[-23,1],pch=19)
text(snum[-23],pve[-23,1],1:23)

