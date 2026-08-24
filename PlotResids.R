source("plotting_funcs.R")

pdf("sfig_resids.pdf",width=9,height=9)
par(mfrow=c(5,5))
par(mar=c(2,3,2,1))
m<-0
for(j in 1:23){

    prefix<-paste("tro_ch",j,"_m",m,sep="")
    R<-plot_resid2(prefix,"poplist.txt",cex=.3)
    title(main=j,cex=.65)
}
dev.off()



