## Map and PCA (tree added seperately)

## for a pretty map
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)

library(data.table)

## get the allele frequency data

a1f<-list.files(pattern="ad1_fff")
a2f<-a1f
a2f<-gsub("ad1","ad2",a2f)
N<-length(a1f)
ids<-read.table("IDs.txt",header=FALSE)
temp<-gsub("ad1_fff_lycSpecPool_chrom","",a1f)
chrom<-gsub(".txt","",temp)
reps<-grep(pattern="rep",ids[,1])

P<-vector("list",24)
n<-vector("list",24)
for(i in 1:N){
	a1<-as.matrix(fread(a1f[i],header=F))
	a2<-as.matrix(fread(a2f[i],header=F))
	n[[i]]<-a1+a2
	P[[i]]<-a2/(a1+a2) ## non-ref
}

## generate the map
locations<-read.table("mapIds.txt",header=TRUE)
# Convert to spatial object
locations_sf <- st_as_sf(locations, coords = c("lon", "lat"), crs = 4326)

world <- ne_countries(scale = "medium", returnclass = "sf")
states <- ne_states(country = "united states of america", returnclass = "sf")

western_usa_map <- ggplot() +
  geom_sf(data = world %>% filter(admin == "United States of America"), fill = "gray95", color = "gray60") +
  geom_sf(data = states, color = "gray80", fill = NA) +
  geom_sf(data = locations_sf, aes(color = color), size = 3, show.legend = FALSE) +
  theme_minimal() +
  coord_sf(xlim = c(-130, -105), ylim = c(30, 50))# +
 # labs(title = "(A) Map")

pdf("fig_mapPCA.pdf",width=8,heigh=8)
cm<-1.4;ca<-1.1;cl<-1.4
par(mfrow=c(2,2))
par(mar=c(4,5,2.5,1.5))
plot(1,1,type='n',axes=FALSE,xlab="",ylab="")
title(main="(a) Locality map",cex.main=cm)
plot(1,1,type='n',axes=FALSE,xlab="",ylab="")
title(main="(b) Chronogram",cex.main=cm)

## PC auto, exluces Z and 21.2
Pa<-as.matrix(rbind(P[[1]],P[[2]],P[[3]],P[[4]],P[[5]],P[[6]],P[[7]],P[[8]],P[[9]],
		    P[[10]],P[[1]],P[[12]],P[[13]],P[[14]],P[[16]],P[[18]],P[[19]],
		    P[[20]],P[[21]],P[[22]],P[[23]],P[[24]]))
Pa[is.na(Pa)]<-0.01
pc<-prcomp(t(Pa[,-c(reps,17)]),center=TRUE,scale=FALSE) ## drop replicates and MEN12
o<-summary(pc)
pct<-round(o$importance[2,1:3] * 100,1)

plot(-1*pc$x[,1],-1*pc$x[,2],pch=19,col=ids[-c(reps,17),2],xlab=paste("PC1 (",pct[1],")",sep=""),ylab=paste("PC2 (",pct[2],")",sep=""),cex.lab=cl,cex.axis=ca)
title(main="(c) PCA autosomes",cex.main=cm)
 
pops<-gsub(pattern="[0-9]+",replacement="",gsub(pattern="Lyc-",replacement="",ids[,1]))
legend(-400,420,pops[-c(reps,17)],pch=19,col=ids[-c(reps,17),2],ncol=3,cex=.6,bty='n')



## Z
Pz<-P[[17]]
Pz[is.na(Pz)]<-0.01
pcz<-prcomp(t(Pz[,-c(reps,17)]),center=TRUE,scale=FALSE) ## drop replicates and MEN12
o<-summary(pcz)
pct<-round(o$importance[2,1:3] * 100,1)

plot(-1*pcz$x[,1],pcz$x[,2],pch=19,col=ids[-c(reps,17),2],xlab=paste("PC1 (",pct[1],")",sep=""),ylab=paste("PC2 (",pct[2],")",sep=""),cex.lab=cl,cex.axis=ca)
title(main="(d) PCA Z chromosome",cex.main=cm)

dev.off()
