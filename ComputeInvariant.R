## determines the number of invariant bases to add to the beast XML file

## read in base counts (A, C, G, T) by scaffold from the genome 
## get this with: perl countBases.pl > baseCounts.txt, from /uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome/

cnts<-read.table("baseCounts.txt",header=FALSE)
## get big scaffolds, top 23
totals<-apply(cnts[,-1],1,sum)
rev(sort(totals))[1:23]
# [1] 31424484 26391115 25798047 25285067 25079126 25067872 22431864 21864934
# [9] 21544644 21226148 20022703 19681008 18787320 18503979 18445203 17279337
#[17] 17086116 16627315 16626411 16602730 15714686 13063727  9211676

chr<-which(totals >= 9211676)

## get total A, C, G, T
bcnt<-apply(cnts[chr,-1],2,sum)

## multiply by 0.00035 to match subsetting for SNP data
prop<-0.00035
sbcnt<-floor(bcnt*prop)

## get SNP bases
## perl countBases.pl > snpCounts.txt
snps<-read.table("snpCounts.txt",header=FALSE)
## average across Lycaeides
snpCnts<-floor(apply(snps[-13,-1],2,mean))
invar<-sbcnt-snpCnts
invar
#   A    C    G    T 
#50190 28322 28284 50114

