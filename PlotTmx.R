library(scales)
library(rstan)
library(loo)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## load summaries from summariz treemix
load(file="tmx.rdat")

source("plotting_funcs.R")

plot_tree_internal <- function(d, e, o = NA, cex = 1, disp = 0.005, plus = 0.005,
                               arrow = 0.05, ybar = 0.01, scale = FALSE, mbar = TRUE,
                               mse = 0.01, plotmig = TRUE, plotnames = TRUE,
                               xmin = 0, lwd = 1, font = 1) {
  plot(d$x, d$y, axes = FALSE, ylab = "", xlab = "Drift parameter",
       xlim = c(xmin, max(d$x) + plus), pch = "",cex.lab=1.4)
  axis(1)

  maxx<-max(d$x) + plus

  mw <- max(e[e[, 5] == "MIG", 4], na.rm = TRUE)
 # mcols <- rev(heat.colors(150))
  mcols <- c(rep(NA,49),alpha("firebrick",seq(0,1,length.out=101)))



  for (i in seq_len(nrow(e))) {
    col <- "black"
    if (e[i, 5] == "MIG") {
      w <- floor(e[i, 4] * 200) + 50
      if (mw > 0.5) w <- floor(e[i, 4] * 100) + 50
      col <- mcols[w]
      if (is.na(col)) col <- "blue"
    }

    v1 <- d[d[, 1] == e[i, 1], ]
    v2 <- d[d[, 1] == e[i, 2], ]

    if (e[i, 5] == "MIG") {
      if (plotmig) {
	#arrows( v1[1,]$x, v1[1,]$y, v2[1,]$x, v2[1,]$y, col = col, length = arrow)
        arrows(v1$x, v1$y, v2$x, v2$y, col = col, length = arrow)
      }
    } else {
	lines( c(v1[1,]$x, v2[1,]$x), c(v1[1,]$y, v2[1,]$y), col = col, lwd = lwd)
    
     ## tried to get orthognoal branches, works but messes up the precise placement of
     ## the migration edges  
     # Orthogonal branches: vertical then horizontal
     # x0 <- v1$x
     # y0 <- v1$y
     # x1 <- v2$x
     # y1 <- v2$y
     # lines(c(x0, x0), c(y0, y1), col = col, lwd = lwd)
     # lines(c(x0, x1), c(y1, y1), col = col, lwd = lwd)
    }
  }

  tmp <- d[d[, 5] == "TIP", ]
  if (!is.na(o)[1]) {
    for (i in seq_len(nrow(tmp))) {
      tcol <- o[o[, 1] == tmp[i, 2], 2]
      if (plotnames) {
        text(tmp[i, ]$x + disp, tmp[i, ]$y, labels = tmp[i, 2], adj = 0,
             cex = cex, col = tcol, font = font)
      }
    }
  } else if (plotnames) {
    tlabs<-gsub(pattern="[0-9]+",replacement="",x=tmp[,2])
    text(rep(maxx-.75*disp,length(tmp$x)), tmp$y, labels = tlabs, adj = 0, cex = cex*.9, font = font)
    for(i in 1:length(tmp$x)){ 
	  lines(c(tmp$x[i],maxx-disp), c(tmp$y[i],tmp$y[i]), lty=3)
    }
  }

  if (scale) {
    lines(c(0, mse * 10), c(ybar, ybar))
    text(0, ybar - 0.04, lab = "10 s.e.", adj = 0, cex = 0.8)
    lines(c(0, 0), c(ybar - 0.01, ybar + 0.01))
    lines(c(mse * 10, mse * 10), c(ybar - 0.01, ybar + 0.01))
  }

  if (mbar) {
    mcols <- alpha("firebrick",seq(0,1,length.out=101))
    #mcols <- rev(heat.colors(150))[50:150]
    ymi <- ybar + 0.15
    yma <- ybar + 0.35
    l <- 0.2
    w <- l / 100
    xma <- max(d$x / 20)
    rect(rep(xma+disp, 100), ymi + (0:99) * w, rep(xma, 100), ymi + (1:100) * w,
         col = mcols, border = mcols)
    text(xma + disp*1.4, ymi, lab = "0", adj = 0, cex = 0.7)
    if (mw > 0.5) {
      text(xma + disp*1.4, yma, lab = "1", adj = 0, cex = 0.7)
    } else {
      text(xma + disp*1.4, yma, lab = "0.5", adj = 0, cex = 0.7)
    }
    text(xma, yma + 0.09, lab = "Admixture", adj = 0, cex = 0.85)
    text(xma, yma + 0.04, lab = "proportion", adj = 0, cex = 0.85)
  }
}

plot_treemix_tree <- function(stem, mbar = TRUE) {
  # Read TreeMix output files
  vertices <- read.table(gzfile(paste0(stem, ".vertices.gz")), as.is = TRUE)
  edges <- read.table(gzfile(paste0(stem, ".edges.gz")), as.is = TRUE)
  covse <- read.table(gzfile(paste0(stem, ".covse.gz")), as.is = TRUE)

  edges[,3] = edges[,3]*edges[,4]
  edges[,3] = edges[,3]*edges[,4]


  # Estimate average SE for scale bar
  mse <- mean(apply(covse, 1, mean))

  # Clean and prep vertices
  vertices$x <- NA
  vertices$y <- NA
  vertices$ymin <- NA
  vertices$ymax <- NA

  # Apply TreeMix logic
  vertices <- set_y_coords(vertices)
  vertices <- set_x_coords(vertices, edges)
  vertices <- set_mig_coords(vertices, edges)

  # Plot
  plot_tree_internal(vertices, edges, mse = mse, mbar = mbar)
}

## example plots
pdf("tp_ch2_m2.pdf",width=5,height=5)
par(mar=c(4,1,1,2.5))
plot_treemix_tree("tro_ch2_m2")
dev.off()


## has scaffold, chrom#, chrom size (bp), number snps
ssz<-read.table("ChromNameSize.txt",header=FALSE)

## gene density, also includes chrom sizes
gd<-read.table("GeneDensitySize.txt",header=TRUE)


rho<-read.table("../Recomb_iSMC/ChromRho.txt",header=TRUE)

logpve<-log(pve/(1-pve))## deals with the asymtote towards 100% pve
m<-seq(0,8,1)
betas<-rep(NA,22)
for(i in 1:22){
	o<-lm(logpve[i,] ~ m)
	betas[i]<-o$coefficients[2]
}	

################ matrix version ###########

## loop over models
## overall design matrix
stand<-function(x){
	xx<-(x-mean(x))/sd(x)
	return(xx)
}

Xmain<-as.matrix(cbind(rep(1,22),stand(ssz[-23,3]),stand(rho$mnRho),stand(gd$gden[-23])))

pc<-prcomp(Xmain[,-1],center=TRUE,scale=TRUE)
pc
#Standard deviations (1, .., p=3):
#[1] 1.4346269 0.8359448 0.4929928

#Rotation (n x k) = (3 x 3):
#            PC1        PC2         PC3
#[1,] -0.6219061 -0.3471719 -0.70192914
#[2,]  0.6257294  0.3186393 -0.71199140
#[3,] -0.4708456  0.8820095 -0.01907195
## PC1, bigger values = shorter chromosomes, with more recombination and fewer genes
## flip it
summary(pc)
#Importance of components:
#                          PC1    PC2     PC3
#Standard deviation     1.4346 0.8359 0.49299
#Proportion of Variance 0.6861 0.2329 0.08101
#Cumulative Proportion  0.6861 0.9190 1.0000

Xpc<-as.matrix(cbind(rep(1,22),-1*stand(pc$x[,1])))

##########################################

sfunc<-function(Y=y,X=NULL,int=FALSE){
	if(is.null(X)==TRUE){
		X<-as.matrix(rep(1,length(Y)))
		form<-Y~.
	}else if(int==FALSE){
		form<-Y~.
	} else{
		form<-Y~.^2
	}
	df<-data.frame(Y=Y,X=X)
	X<-model.matrix(form,data=df)
	dat<-list(y=Y,x=as.matrix(X),K=dim(X)[2],N=dim(X)[1])
	o<-stan("LM.stan",data=dat,seed=31722)
	ll_matrix<-extract_log_lik(o)
	ol<-loo(ll_matrix)
	return(ol$estimates[1,])
}

dpve<-dim(pve)[2] ## 9
elpd<-matrix(NA,nrow=13,ncol=9)
elpdse<-matrix(NA,nrow=13,ncol=9)
des<-matrix(c(2,3,2,4,3,4),nrow=3,byrow=TRUE)

for(i in 1:dpve){
	ypve<-stand(pve[-23,i])
	## null
	out<-sfunc(Y=ypve)
	elpd[1,i]<-out[1]
	elpdse[1,i]<-out[2]
	## singles
	for(j in 2:4){
		out<-sfunc(Y=ypve,X=Xmain[,j])
		elpd[j,i]<-out[1]
		elpdse[j,i]<-out[2]
	}
	## two ways, no int
	for(j in 1:3){
		out<-sfunc(Y=ypve,X=Xmain[,des[j,]])
		elpd[4+j,i]<-out[1]
		elpdse[4+j,i]<-out[2]
	}
	## two ways, int
	for(j in 1:3){
		out<-sfunc(Y=ypve,X=Xmain[,des[j,]],int=TRUE)
		elpd[7+j,i]<-out[1]
		elpdse[7+j,i]<-out[2]
	}
	## three way, and three way int
	out<-sfunc(Y=ypve,X=Xmain[,2:4])
	elpd[11,i]<-out[1]
	elpdse[11,i]<-out[2]
	out<-sfunc(Y=ypve,X=Xmain[,2:4],int=TRUE)
	elpd[12,i]<-out[1]
	elpdse[12,i]<-out[2]
	## PC version
	out<-sfunc(Y=ypve,X=Xpc[,2],int=TRUE)
	elpd[13,i]<-out[1]
	elpdse[13,i]<-out[2]
}

mods<-c("null","Sz","Rh","Gd","SzRh","SzGd","RhGd","SzRhI","SzGdI","RhGdI","SzRhGd","SzRhGdI")
## best model by m = 0 .. 9
mods[apply(elpd[1:12,],2,which.max)]
#[1] "Gd"    "SzRhI" "null"  "null"  "SzRh"  "Sz"    "Sz"    "Sz"    "Rh" 

## best model across all migration rates
mods[which.max(apply(elpd[1:12,],1,sum))]
#[1] "Sz"

## with pc
pcrnk<-rep(NA,9)
for(i in 1:9){
	pcrnk[i]<-sum(elpd[1:12,i] > elpd[13,i]) + 1
}
pcrnk
#[1] 2 7 6 3 7 2 2 1 2
median(pcrnk)
#[1] 2
mean(pcrnk)
#[1] 3.555556

## stan fit
so<-vector("list",9)
for(i in 1:9){
        ypve<-stand(pve[-23,i])
	dat<-list(y=ypve,x=Xpc,K=2,N=dim(Xpc)[1])
	so[[i]]<-stan("LM.stan",data=dat,seed=31722)
}

pdf("sfig_tmixpve.pdf",width=9,height=9)
cl<-1.4;cm<-1.4
par(mfrow=c(3,3))
par(mar=c(4.5,4.5,2.5,1))
for(i in 1:9){
        ypve<-stand(pve[-23,i])
	dat<-list(y=ypve,x=Xpc,K=2,N=dim(Xpc)[1])
	parms<-extract(so[[i]])
	## get posterior for plot
	xx<-seq(-3,3,0.01)
	yhat<-matrix(nrow=4000,ncol=length(xx)) ## 4000 HMC samples
	for(k in 1:4000){
	        yhat[k,]<-parms[[1]][k,1] + parms[[1]][k,2]*xx
	}
	est<-apply(yhat,2,quantile,probs=c(.5,.05,.95))
	plot(dat$x[,2],dat$y,axes=FALSE,xlab="Chromosome PC1",ylab="Graph PCvE",cex.lab=cl,type='n')
	polygon(c(xx,rev(xx)),c(est[2,],rev(est[3,])),col=alpha("violet",.3),border=NA)
	points(dat$x[,2],dat$y,pch=19)
	lines(xx,est[1,],col="violet",lwd=1.5)
	ya<-seq(round(min(pve[-23,i]),2),round(max(pve[-23,i]),3),length.out=4)
	axis(2,at=(ya-mean(pve[-23,i]))/sd(pve[-23,i]),round(ya,3))
	axis(1)
	box()
	title(main=paste("(",LETTERS[i],") ",i," migration edges",sep=""),cex.main=cm)
}

dev.off()

est<-matrix(NA,nrow=9,ncol=5)
for(i in 1:9){
	parm<-extract(so[[i]])
	est[i,]<-quantile(parm[[1]][,2],probs=c(.5,.05,.95,.25,.75))
}

cor.test(0:8,est[,1])
#	Pearson's product-moment correlation
#
#data:  0:8 and est[, 1]
#t = 5.1948, df = 7, p-value = 0.00126
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.5559971 0.9770183
3sample estimates:
#      cor 
#0.8910856 


pdf("sfig_cors.pdf",width=6,height=4.5)
dev.off()

save(list=ls(),file="tmixAlt.rdat")

pdf("fig_tmix_alt.pdf",width=9,height=9)
cl<-1.4;cm<-1.4;ca<-1
layout(matrix(c(1,2,2,3,4,5,5,6,7,8,8,9),nrow=3,ncol=4,byrow=TRUE),widths=c(3,1.5,1.5,3))
par(mar=c(4.5,4.5,2.5,1))

plot(seq(0,8,1),pve[1,],ylim=c(.83,1),type='l',xlab="Admixture edges",ylab="Graph PCvE",col="gray60",cex.lab=cl)
for(i in 2:22){lines(seq(0,8,1),pve[i,],col="gray60")}
lines(seq(0,8,1),pve[23,],lwd=2)
legend(5,.885,c("autosome","Z chromosome"),lwd=c(1,1.5),col=c("gray60","black"),bty='n')
title(main="(A) Allele frequency covariance explained",cex.main=cm)

plot(pve[,1],eigenRank[,1],pch=19,xlab="Tree PCvE",ylab="Eigen Rank",cex.lab=cl,cex.axis=ca,col=c(rep("gray60",22),"black"))
title(main="(B) Tree variance explained and residual structure",cex.main=cm)

plot(0:8,est[,1],pch=19,ylim=c(-1,1),xlab="Admixture edges",ylab="Coefficient",cex.lab=1.3,cex.axis=1.1)
segments(0:8,est[,4],0:8,est[,5],lwd=4.5)
segments(0:8,est[,2],0:8,est[,3])
abline(h=0,lty=2)
title(main="(C) Graph PcVE and chromosome features",cex.main=cm)


dev.off() ## do not need the rest

####### original version 2merge results in figure#####

## center and standardize covariate and response
ypve<-(pve[-23,1]-mean(pve[-23,1]))/sd(pve[-23,3])
xssz<-(ssz[-23,3]-mean(ssz[-23,3]))/sd(ssz[-23,3])
dat<-list(y=ypve,x=xssz,N=length(xssz))
o<-stan("lm.stan",data=dat)
## get posterior
parms<-extract(o)
## get posterior for plot
xx<-seq(-3,3,0.01)
yhat<-matrix(nrow=4000,ncol=length(xx)) ## 4000 HMC samples
for(i in 1:4000){
	yhat[i,]<-parms[[1]][i] + parms[[2]][i]*xx
}
est<-apply(yhat,2,quantile,probs=c(.5,.05,.95))

##############################################

pdf("fig_tmix.pdf",width=9,height=9)
cl<-1.4;cm<-1.4
layout(matrix(c(1,1,2,2,3,4,4,5,6,7,7,8),nrow=3,ncol=4,byrow=TRUE),widths=c(3,1.5,1.5,3))
par(mar=c(4.5,4.5,2.5,1))

plot(seq(0,8,1),pve[1,],ylim=c(.83,1),type='l',xlab="Admixture edges",ylab="Graph PCVE",col="gray60",cex.lab=cl)
for(i in 2:22){lines(seq(0,8,1),pve[i,],col="gray60")}
lines(seq(0,8,1),pve[23,],lwd=2)
legend(5,.885,c("autosome","Z chromosome"),lwd=c(1,1.5),col=c("gray60","black"),bty='n')
title(main="(A) Allele frequency covariance explained",cex.main=cm)

plot(dat$x,dat$y,axes=FALSE,xlab="Chromsome length (Mbps)",ylab="Graph PCVE",cex.lab=cl)
polygon(c(xx,rev(xx)),c(est[2,],rev(est[3,])),col=alpha("violet",.3),border=NA)
points(dat$x,dat$y,pch=19)
lines(xx,est[1,],col="violet",lwd=1.5)
ya<-seq(round(min(pve[-23,1]),2),round(max(pve[-23,1]),2),by=0.01)
xa<-seq(from=10,to=33,by=3)
axis(2,at=(ya-mean(pve[-23,1]))/sd(pve[-23,1]),ya)
axis(1,at=(xa-mean(ssz[-23,3]/1000000))/sd(ssz[-23,3]/1000000),xa)
box()
olm<-lm(dat$y ~ dat$x)
oo<-summary(olm)
text(-1.8,-2,expression(paste(r^2," = 0.12")),cex=1.25)
#text(-2,-2,expression(paste(r^2," = ",round(oo$r.squared,2))))
title(main="(B) Chromsome size vs PCVE",cex.main=cm)

par(mar=c(4,1,2.5,2.5))
plot_treemix_tree("tro_ch2_m0",mbar=FALSE)
title(main="(C) Chromosome 2 m=0",cex.main=cm)
plot_treemix_tree("tro_ch4_m0",mbar=FALSE)
title(main="(D) Chromosome 4 m=0",cex.main=cm)
plot_treemix_tree("tro_ch23_m0",mbar=FALSE)
title(main="(E) Chromosome Z m=0",cex.main=cm)
plot_treemix_tree("tro_ch2_m2",mbar=FALSE)
title(main="(F) Chromosome 2 m=2",cex.main=cm)
plot_treemix_tree("tro_ch4_m2",mbar=FALSE)
title(main="(G) Chromosome 4 m=2",cex.main=cm)
plot_treemix_tree("tro_ch23_m2",mbar=TRUE)
title(main="(H) Chromosome Z m=2",cex.main=cm)
dev.off()

## using rho instead ################3
ypve<-(pve[-23,1]-mean(pve[-23,1]))/sd(pve[-23,3])

k<-2
xrho<-(rho[,k]-mean(rho[,k]))/sd(rho[,k])
summary(lm(ypve~xrho))

dat<-list(y=ypve,x=xrho,N=length(xrho))
o<-stan("lm.stan",data=dat)
## get posterior
parms<-extract(o)
## get posterior for plot
xx<-seq(-1.8,1.8,0.01)
yhat<-matrix(nrow=4000,ncol=length(xx)) ## 4000 HMC samples
for(i in 1:4000){
	yhat[i,]<-parms[[1]][i] + parms[[2]][i]*xx
}
est<-apply(yhat,2,quantile,probs=c(.5,.05,.95))

pdf("fig_tmix_rho.pdf",width=9,height=9)
cl<-1.4;cm<-1.4
layout(matrix(c(1,1,2,2,3,4,4,5,6,7,7,8),nrow=3,ncol=4,byrow=TRUE),widths=c(3,1.5,1.5,3))
par(mar=c(4.5,4.5,2.5,1))

plot(seq(0,8,1),pve[1,],ylim=c(.83,1),type='l',xlab="Admixture edges",ylab="Graph PCVE",col="gray60",cex.lab=cl)
for(i in 2:22){lines(seq(0,8,1),pve[i,],col="gray60")}
lines(seq(0,8,1),pve[23,],lwd=2)
legend(5,.885,c("autosome","Z chromosome"),lwd=c(1,1.5),col=c("gray60","black"),bty='n')
title(main="(A) Allele frequency covariance explained",cex.main=cm)

plot(dat$x,dat$y,axes=FALSE,xlab="Recombination rate (rho)",ylab="Graph PCVE",cex.lab=cl)
polygon(c(xx,rev(xx)),c(est[2,],rev(est[3,])),col=alpha("violet",.3),border=NA)
points(dat$x,dat$y,pch=19)
lines(xx,est[1,],col="violet",lwd=1.5)
ya<-seq(round(min(pve[-23,1]),2),round(max(pve[-23,1]),2),by=0.01)
xa<-seq(0.123,0.141,length.out=10)
axis(2,at=(ya-mean(pve[-23,1]))/sd(pve[-23,1]),ya)
axis(1,at=(xa-mean(rho[,2]))/sd(rho[,2]),xa)
box()
olm<-lm(dat$y ~ dat$x)
oo<-summary(olm)
text(-1.5,-2,expression(paste(r^2," = 0.05")),cex=1.25)
#text(-2,-2,expression(paste(r^2," = ",round(oo$r.squared,2))))
title(main="(B) Recombination vs PCVE",cex.main=cm)

par(mar=c(4,1,2.5,2.5))
plot_treemix_tree("tro_ch2_m0",mbar=FALSE)
title(main="(C) Chromosome 2 m=0",cex.main=cm)
plot_treemix_tree("tro_ch4_m0",mbar=FALSE)
title(main="(D) Chromosome 4 m=0",cex.main=cm)
plot_treemix_tree("tro_ch23_m0",mbar=FALSE)
title(main="(E) Chromosome Z m=0",cex.main=cm)
plot_treemix_tree("tro_ch2_m2",mbar=FALSE)
title(main="(F) Chromosome 2 m=2",cex.main=cm)
plot_treemix_tree("tro_ch4_m2",mbar=FALSE)
title(main="(G) Chromosome 4 m=2",cex.main=cm)
plot_treemix_tree("tro_ch23_m2",mbar=TRUE)
title(main="(H) Chromosome Z m=2",cex.main=cm)
dev.off()

## si fig, all the graphs
for(i in 1:23){
	of<-paste("sfig_tmx_ch",i,".pdf",sep="")
	pdf(of,width=9,height=6)
	par(mfrow=c(2,3))
	par(mar=c(4.5,4.5,2.5,1))
	for(j in seq(0,5,1)){
		ff<-paste("tro_ch",i,"_m",j,sep="")
		if(j==5){
			plot_treemix_tree(ff,mbar=TRUE)
		} else{
			plot_treemix_tree(ff,mbar=FALSE)
		}
	}
	dev.off()
}

## rest of regressions
## m = 0
#Inference for Stan model: anon_model.
#4 chains, each with iter=2000; warmup=1000; thin=1;
#post-warmup draws per chain=1000, total post-warmup draws=4000.

#        mean se_mean   sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
#beta0   0.00    0.00 0.26  -0.51  -0.18   0.01   0.18   0.52  3655    1
#beta1  -0.41    0.00 0.26  -0.90  -0.59  -0.42  -0.24   0.11  2809    1
#sigma   1.21    0.00 0.21   0.88   1.06   1.19   1.33   1.69  2762    1
#lp__  -14.58    0.03 1.29 -17.86 -15.18 -14.25 -13.65 -13.09  1867    1
obm<-vector("list",9)
obm[[1]]<-o
for(i in 2:9){
	ypve<-(pve[-23,i]-mean(pve[-23,i]))/sd(pve[-23,i])
	dat<-list(y=ypve,x=xssz,N=length(xssz))
	o<-stan("lm.stan",data=dat)
	obm[[i]]<-o
}

est<-matrix(NA,nrow=9,ncol=5)
for(i in 1:9){
	parm<-extract(obm[[i]])
	est[i,]<-quantile(parm$beta1,probs=c(.5,.025,.975,.25,.75))
}

cor.test(0:8,est[,1])
#Pearson's product-moment correlation

#data:  0:8 and est[, 1]
#t = 3.519, df = 7, p-value = 0.009741
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.2880982 0.9559622
#sample estimates:
#    cor
#0.79929

pdf("sfig_cors.pdf",width=6,height=4.5)
par(mar=c(5,5,1,1))
plot(0:8,est[,1],pch=19,ylim=c(-1,1),xlab="Migration edges",ylab="Coefficient",cex.lab=1.3,cex.axis=1.1)
segments(0:8,est[,4],0:8,est[,5],lwd=4.5)
segments(0:8,est[,2],0:8,est[,3])
abline(h=0,lty=2)
dev.off()

## rest of regressions, with rho
obm<-vector("list",9)
obm[[1]]<-o
for(i in 2:9){
	ypve<-(pve[-23,i]-mean(pve[-23,i]))/sd(pve[-23,i])
	dat<-list(y=ypve,x=xrho,N=length(xssz))
	o<-stan("lm.stan",data=dat)
	obm[[i]]<-o
}

est<-matrix(NA,nrow=9,ncol=5)
for(i in 1:9){
	parm<-extract(obm[[i]])
	est[i,]<-quantile(parm$beta1,probs=c(.5,.05,.95,.25,.75))
}

cor.test(0:8,est[,1])
#Pearson's product-moment correlation
#	Pearson's product-moment correlation
#
#data:  0:8 and est[, 1]
#t = -6.3023, df = 7, p-value = 0.0004033
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.9837610 -0.6653567
#sample estimates:
#       cor 
#-0.9220465 


pdf("sfig_cors_rho.pdf",width=6,height=4.5)
par(mar=c(5,5,1,1))
plot(0:8,est[,1],pch=19,ylim=c(-.8,.8),xlab="Migration edges",ylab="Coefficient",cex.lab=1.3,cex.axis=1.1)
segments(0:8,est[,4],0:8,est[,5],lwd=4.5)
segments(0:8,est[,2],0:8,est[,3])
abline(h=0,lty=2)
dev.off()

## get migration edge weights
edgeweights<-function(stem) {
  # Read TreeMix output files
  edges <- read.table(gzfile(paste0(stem, ".edges.gz")), as.is = TRUE)
  return(edges$V4[edges$V5=="MIG"])
}

j<-1
admix1<-matrix(NA,nrow=23,ncol=1)
for(i in 1:23){
	ff<-paste("tro_ch",i,"_m",j,sep="")
	admix1[i,]<-edgeweights(ff)
}
summary(admix1[-23,])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2696  0.3557  0.3948  0.3896  0.4218  0.4795 
#[23,] 0.436630
j<-2
admix2<-matrix(NA,nrow=23,ncol=2)
for(i in 1:23){
	ff<-paste("tro_ch",i,"_m",j,sep="")
	admix2[i,]<-edgeweights(ff)
}
summary(as.vector(admix2[-23,]))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.16edgeweights(ff)
13  0.3161  0.3751  0.3687  0.4264  0.4970 
#[23,] 0.296775 0.362413
j<-6
admix6<-matrix(NA,nrow=23,ncol=6)
for(i in 1:23){
	ff<-paste("tro_ch",i,"_m",j,sep="")
	admix6[i,]<-edgeweights(ff)
}
summary(as.vector(admix6[-23,]))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.05834 0.25336 0.34497 0.33588 0.41742 0.61240 
admix6[23,]
#[1] 0.294562 0.455030 0.180839 0.428366 0.134007 0.199130


save(list=ls(),file="tmxFigs.rdat")
