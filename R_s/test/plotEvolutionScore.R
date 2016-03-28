setwd("/home/limeng/splicingSNP/R/")
library(ggplot2)
source("multiPlot.r")

setwd("/home/limeng/splicingSNP/code/")


dataAllExons<-read.table("allExons.rand.csv.phylop",header=FALSE,as.is=TRUE);
dataAllExons<-cbind(dataAllExons[,1],rep("all.exons",nrow(dataAllExons)));
dataAllExons<-unname(dataAllExons);


#splicingJunction.20.highPSI<-read.table("high.psi.exons.csv.phylop",as.is=TRUE,header=FALSE);
splicingJunction.20.highPSI<-read.table("10high.psi.exons.csv.phylop",as.is=TRUE,header=FALSE);
splicingJunction.20.highPSI<-cbind(splicingJunction.20.highPSI[,1],rep("splicingJunction.20.highPSI",nrow(splicingJunction.20.highPSI)));
splicingJunction.20.highPSI<-unname(splicingJunction.20.highPSI)


#splicingJunction.20.lowPSI<-read.table("low.psi.exons.csv.phylop",as.is=TRUE,header=FALSE)
splicingJunction.20.lowPSI<-read.table("10low.psi.exons.csv.phylop",as.is=TRUE,header=FALSE)
splicingJunction.20.lowPSI<-cbind(splicingJunction.20.lowPSI[,1],rep("splicingJunction.20.lowPSI",nrow(splicingJunction.20.lowPSI)));
splicingJunction.20.lowPSI<-unname(splicingJunction.20.lowPSI)


hgmd<-read.table("hgmdNearestExonRegions25.unique.csv.phylop",header=FALSE,as.is=TRUE)
#hgmd<-read.table("hgmdNearestExonRegions25.csv.phylop",header=FALSE,as.is=TRUE)
#hgmd<-read.table("hgmdNearestExonRegionsNonunique.csv.phylop",header=FALSE,as.is=TRUE);
hgmd<-cbind(hgmd[,1],rep("hgmd",nrow(hgmd) )  );
hgmd<-unname(hgmd);


#splicingJunction.20.AF10<-read.table("splicingJunction20.AF10.exonOnly.exonNearSNP.phylop")
#splicingJunction.20.AF10<-cbind(splicingJunction.20.AF10[,1],rep("splicingJunction.20.AF10",nrow(splicingJunction.20.AF10)))
#splicingJunction.20.AF10<-unname(splicingJunction.20.AF10)


allData<-rbind(dataAllExons,splicingJunction.20.highPSI,hgmd,splicingJunction.20.lowPSI);
allData<-data.frame(evolution=as.numeric(allData[,1]),label=allData[,2]) ;
colnames(allData)<-c("evolution.score","label")


cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# To use for line and point colors, add
#scale_colour_manual(values=cbPalette)


pd<-ggplot(allData,aes(x=evolution.score,color=label))
pd<-pd+geom_density()
pd<-pd+scale_colour_manual(values=cbPalette)
pd<-pd+xlim(-0.1,0.6)
#+geom_histogram(position="dodge",aes(y = ..density..))
#print(p)


ph<-ggplot(allData,aes(x=evolution.score,fill=label))
ph<-ph+scale_colour_manual(values=cbPalette)
ph<-ph+xlim(-0.1,0.6)
ph<-ph+geom_histogram(position="dodge",aes(y = ..density..))
#print(ph)


phHighPSIAndHGMD<-ggplot(data=subset(allData,(label=="splicingJunction.20.highPSI")|(label=="hgmd")),aes(x=evolution.score,fill=label))
phHighPSIAndHGMD<-phHighPSIAndHGMD+scale_colour_manual(values=cbPalette)
phHighPSIAndHGMD<-phHighPSIAndHGMD+xlim(-0.1,0.6)
phHighPSIAndHGMD<-phHighPSIAndHGMD+geom_histogram(position="dodge",aes(y = ..density..))
print(phHighPSIAndHGMD)
#print(ph)


pdHighPSIAndHGMD<-ggplot(subset(allData,(label=="splicingJunction.20.highPSI")|(label=="hgmd")),aes(x=evolution.score,color=label))
pdHighPSIAndHGMD<-pdHighPSIAndHGMD+scale_colour_manual(values=cbPalette)
pdHighPSIAndHGMD<-pdHighPSIAndHGMD+xlim(-0.1,0.6)
pdHighPSIAndHGMD<-pdHighPSIAndHGMD+geom_density();
#print(ph)


multiplot(pd,ph,pdHighPSIAndHGMD,phHighPSIAndHGMD,cols=2)






















