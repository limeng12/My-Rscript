

mBed=new("BED");

mBed<-readFromBedFile(b,"/home/limeng/splicingSNP/annotation/refgene.bed");

#exonPos<-getExonRegion(mBed,"NM_001918",1);

allExons<-getAllExonCoordinate(mBed);

#onlyUsefullExons<-allExons[grepl("chr[1-9,x,y]:|chr1[0-9]:|chr2[0-2]:",allExons)];
#cat(str_c("the number of exons in chr1-22,x,y:",length(onlyUsefullExons),"\n"));

setwd("/home/limeng/splicingSNP/code")
cat(allExons,file="allExons.csv",sep="\n");

#1  897738     1  899928
snp<-data.frame(chr=c("1","1","1"),pos=c(897325,897738,899928),stringsAsFactors=FALSE );
targetRegions<-getNearestExons(mBed,snp,3);

exonsCenterPosition<-getNearestExonCenterPosition(mBed,snp,3);

write.csv(convertToDataframe(mBed));






