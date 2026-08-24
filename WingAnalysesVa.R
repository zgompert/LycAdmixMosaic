chd<-read.table("../ChromNameSize.txt",header=FALSE)
#chd<-read.table("../../../../TreeMx/ChromNameSize.txt",header=FALSE)

## compute residual matrix
compRM<-function(x=NA,Y=NA){
	N<-dim(Y)[1]
	R<-Y
	for(j in 1:N){
		o<-lm(Y[j,]~x)
		R[j,]<-o$residuals
	}
	return(R)
}

## test if residuals for Z (23) > null
RMtest<-function(R=NA,nperm=1000,nt=17){
	z<-23
	obs<-mean(R[,z])
	pmat<-matrix(NA,nrow=nt,ncol=2)
	pmat[,1]<-1:nt
	null<-rep(NA,nperm)
	for(k in 1:nperm){
		pmat[,2]<-sample(1:23,nt,replace=FALSE)
		nR<-R[pmat]
		null[k]<-mean(nR)
	}
	out<-c(obs,mean(null),mean(null >= obs))
	return(out)
}

pGNP<-read.table("../gemmaP_GNP.txt",header=FALSE)
#pGNP<-read.table("../../geno/gemmaP_GNP.txt",header=FALSE)

vaGNP<-matrix(NA,nrow=17,ncol=23)
for(i in 1:17){
    ff<-paste("mav_o_lyc_GNP_ph",i,"_ch0.param.txt",sep="")
    mav<-read.table(ff,header=FALSE)
    va<-2*pGNP[,1]*(1-pGNP[,1])*mav[,2]^2
    chr<-as.numeric(regmatches(mav[,1], regexpr("^[0-9]+", mav[,1])))
    vaGNP[i,]<-tapply(X=va,INDEX=chr,sum)
    vaGNP[i,]<-vaGNP[i,]/sum(vaGNP[i,])
}
nSnpsGNP<-table(chr)
YnSnpGNP<-compRM(x=as.numeric(nSnpsGNP),Y=vaGNP)
RMtest(R=YnSnpGNP)
#[1]  0.0510422016 -0.0006891057  0.0240000000

YchrGNP<-compRM(x=chd[,3],Y=vaGNP)
RMtest(R=YchrGNP)
#[1] 0.0289718150 0.0003405563 0.1230000000

RMtest(R=vaGNP)
#[1] 0.08905646 0.04362586 0.04800000

### aurorae only

RMtest(R=YnSnpGNP[1:6,],nt=6)
#[1] 0.088665135 0.001375309 0.086000000

RMtest(R=YchrGNP[1:6,],nt=6)
#[1]  0.0770230527 -0.0000246905  0.0760000000

RMtest(R=vaGNP[1:6,],nt=6)
#[1] 0.13146996 0.04311009 0.08000000

pSIN<-read.table("../gemmaP_SIN.txt",header=FALSE)
#pSIN<-read.table("../../geno/gemmaP_SIN.txt",header=FALSE)

vaSIN<-matrix(NA,nrow=17,ncol=23)
for(i in 1:17){
    ff<-paste("mav_o_lyc_SIN_ph",i,"_ch0.param.txt",sep="")
    mav<-read.table(ff,header=FALSE)
    va<-2*pSIN[,1]*(1-pSIN[,1])*mav[,2]^2
    chr<-as.numeric(regmatches(mav[,1], regexpr("^[0-9]+", mav[,1])))
    vaSIN[i,]<-tapply(X=va,INDEX=chr,sum)
    vaSIN[i,]<-vaSIN[i,]/sum(vaSIN[i,])
}
nSnpsSIN<-table(chr)
YnSnpSIN<-compRM(x=as.numeric(nSnpsSIN),Y=vaSIN)
RMtest(R=YnSnpSIN)
#[1]  0.0051984515 -0.0000153398  0.2050000000

YchrSIN<-compRM(x=chd[,3],Y=vaSIN)
RMtest(R=YchrSIN)
#[1] -1.212303e-02  8.982448e-05  9.270000e-01

RMtest(R=vaSIN)
#[1] 0.04193104 0.04335831 0.37800000

### aurorae only

RMtest(R=YnSnpSIN[1:6,],nt=6)
#[1] -0.0038779583  0.0002835056  0.4970000000
RMtest(R=YchrSIN[1:6,],nt=6)
#[1] -0.018039846 -0.000135967  0.987000000
RMtest(R=vaSIN[1:6,],nt=6)
#[1] 0.03459303 0.04390144 0.71800000

pYBG<-read.table("../gemmaP_YBG.txt",header=FALSE)
#pYBG<-read.table("../../geno/gemmaP_YBG.txt",header=FALSE)

vaYBG<-matrix(NA,nrow=17,ncol=23)
for(i in 1:17){
    ff<-paste("mav_o_lyc_YBG_ph",i,"_ch0.param.txt",sep="")
    mav<-read.table(ff,header=FALSE)
    va<-2*pYBG[,1]*(1-pYBG[,1])*mav[,2]^2
    chr<-as.numeric(regmatches(mav[,1], regexpr("^[0-9]+", mav[,1])))
    vaYBG[i,]<-tapply(X=va,INDEX=chr,sum)
    vaYBG[i,]<-vaYBG[i,]/sum(vaYBG[i,])
}
nSnpsYBG<-table(chr)
YnSnpYBG<-compRM(x=as.numeric(nSnpsYBG),Y=vaYBG)
RMtest(R=YnSnpYBG)
#0.0395644436 0.0004307654 0.0290000000

YchrYBG<-compRM(x=chd[,3],Y=vaYBG)
RMtest(R=YchrYBG)
#[1]  0.0216992048 -0.0004634516  0.1280000000

RMtest(R=vaYBG)
#[1] 0.07486146 0.04284347 0.05500000

### aurorae only

RMtest(R=YnSnpYBG[1:6,],nt=6)
#[1]  0.1170212378 -0.0004239905  0.0100000000

RMtest(R=YchrYBG[1:6,],nt=6)
#[1]  0.0707166586 -0.0001453787  0.0740000000

RMtest(R=vaYBG[1:6,],nt=6)
#[1] 0.14120303 0.04229473 0.05300000
