library(ggplot2)
library(dplyr)
library(plyr)
library(readr)
setwd("/Users/mengli/Documents/projects/splicingSNP/corr/"); 

attr<-read_tsv("Brain_samplesAttributes_1.txt",col_names = F); 
attr_data<-as.data.frame(attr,stringAsFactors=F); 
donors<-substr(attr_data[,2],1,9)
donorNumber<-length(unique(donors));
cat( paste0("donors number: ",donorNumber) ); 

cat( paste0("region number: ",length(unique(attr_data[,14]))) ); 

result<-read.table("region_miso_event_result.tsv",sep="\t",header=T,as.is=T); 
jpeg("sd_fxi_corr.jpeg",width=700,height=500); 
uniqueEventResult<-result[!duplicated(result[,4]),]; 

smoothScatter(uniqueEventResult[,1],uniqueEventResult[,2],
              xlab="standard variation",ylab="Fis",
              main=paste0("cor: ",signif(cor(uniqueEventResult[,2],uniqueEventResult[,1]) ,3) ) ); 
dev.off(); 

jpeg("sd_hist.jpeg"); 
hist(uniqueEventResult[,"psi"],main="distribution of psi",xlab="psi",ylab="frequency"); 
dev.off(); 

jpeg("sd_fxi_bin.jpeg",width=700,height=500); 
labelFis<-sapply(uniqueEventResult[,"fis"],function(x){
  if(x<0.85){return("x<0.85");};
  if(x>=0.85&&x<0.91){return("0.85<=x<0.91")};
  return("x>0.91");
  
});

uniqueEventResult<-cbind(uniqueEventResult,labelFis);
uniqueEventResult[,"labelFis"]<-factor(uniqueEventResult[,"labelFis"]); 
p1<-ggplot(uniqueEventResult)+geom_boxplot(aes(labelFis,sd))+theme_classic();
print(p1);
dev.off();

labelFis<-sapply(result[,"fis"],function(x){
  if(x<0.85){return("fis<0.85");};
  if(x>=0.85&&x<0.91){return("0.85<=fis<0.91")};
  return("fis>0.91");
});
labelPsi<-sapply(result[,"psi"],function(x){
  if(x<=0.2){return("psi<0.2");}
  if(0.2<x&&x<=0.4){return("0.2<psi<=0.4");}
  if(0.4<x&&x<=0.6){return("0.4<psi<=0.6");}
  if(0.6<x&&x<=0.8){return("0.6<psi<=0.8");}
  if(0.8<x&&x<=1){return("0.8<psi<=1")}
}); 
result<-cbind(result,labelFis,labelPsi); 

result[,"labelFis"]<-as.character(result[,"labelFis"]);
result[,"labelPsi"]<-as.character(result[,"labelPsi"]); 

fis_psi<-ddply(result,.(event_name,fis,labelFis),function(x){
  medianPsi<-median(x[,"psi"]);
  return(medianPsi);
  
} ); 

colnames(fis_psi)<-c("event_name","fis","labelFis","median_psi"); 
labelPsi<-sapply(fis_psi[,"median_psi"],function(x){
  if(x<=0.2){return("psi<0.2");}
  if(0.2<x&&x<=0.4){return("0.2<psi<=0.4");}
  if(0.4<x&&x<=0.6){return("0.4<psi<=0.6");}
  if(0.6<x&&x<=0.8){return("0.6<psi<=0.8");}
  if(0.8<x&&x<=1){return("0.8<psi<=1")}
  
});

fis_psi<-cbind(fis_psi,labelPsi); 
fis_psi[,"labelPsi"]<-factor(fis_psi[,"labelPsi"],
                             levels = c("psi<0.2","0.2<psi<=0.4","0.4<psi<=0.6","0.6<psi<=0.8","0.8<psi<=1"),
                             ordered=T ); 


fis_psi[,"labelFis"]<-factor(fis_psi[,"labelFis"],
                             levels=c("fis<0.85","0.85<=fis<0.91","fis>0.91" ),ordered=T); 
jpeg("psi_fis.jpeg",width=900,height=600); 
p3<-ggplot(fis_psi)+geom_density(aes(x=median_psi,color=labelFis))+theme_classic(); 
print(p3); 
dev.off(); 

testP<-wilcox.test(fis_psi[fis_psi[,3]=="fis<0.85","median_psi"],
                   fis_psi[fis_psi[,3]=="0.85<=fis<0.91","median_psi"] )$p.value; 

freByFis<-ddply(fis_psi,.(labelFis),nrow);
fis_psi[,"labelFis"]<-as.character(fis_psi[,"labelFis"]); 
resultPer<-ddply(fis_psi,.(labelPsi),function(x){ 
  fis085Count<-sum(x[,"labelFis"]=="fis<0.85"); 
  fis085091Count<-sum(x[,"labelFis"]=="0.85<=fis<0.91"); 
  fis091Count<-sum(x[,"labelFis"]=="fis>0.91"); 
  allCount<-nrow(x); 
  
  return(c(fis085Count/freByFis[1,"V1"],fis085091Count/freByFis[2,"V1"],
           fis091Count/freByFis[3,"V1"],allCount) ); 
}); 

resultPermelt<-melt(resultPer,c("labelPsi","V4")); 
resultPermelt<-cbind(resultPermelt,resultPermelt[,"V4"]*resultPermelt[,"value"] ); 

colnames(resultPermelt)<-c("labelPsi","allCount","variable","value","count"); 

jpeg("psi_fis_bar.jpeg",width=900,height=600);
p4<-ggplot(resultPermelt,aes(x=labelPsi,y=value,fill=variable) )+
  geom_bar(stat="identity",position='dodge' ,wdith=0.2)+theme_classic()+
  ggtitle(paste0("PSIs' wilcox test between fis<0.85 and 0.85<fis<0.91: P=",signif(testP,3) ))+
  scale_fill_discrete(name="",breaks=c("V1","V2","V3") ,
                      labels=c("Fis<0.85","0.85<Fis<0.91","Fis>0.91")  )+
  ylab("percent")
  #geom_text(aes(label=count),position=position_dodge(width=0.9), );
print(p4)
dev.off(); 

result2<-ddply(result,.(labelPsi,event_name),function(x){ 
  sd_psi<-as.numeric(sd(x[,"psi"]) ); 
  return( c(sd_psi,x[1,"labelFis"],length(sd_psi)) ); 
  
});
colnames(result2)<-c("labelPsi","event_name","sd_psi","labelFis","sample_count"); 

result2[,"sd_psi"]<-as.numeric( result2[,"sd_psi"] ); 
result2<-result2[!is.na(result2[,"sd_psi"]),]; 
result2[,"labelPsi"]<-factor(result2[,"labelPsi"],levels = c("psi<0.2",
"0.2<psi<=0.4","0.4<psi<=0.6","0.6<psi<=0.8","0.8<psi<=1"),ordered=T ); 
result2[,"labelFis"]<-factor(result2[,"labelFis"],
                             levels=c("fis<0.85","0.85<=fis<0.91","fis>0.91" ),ordered=T ); 

eventCountResult<-ddply(result2, .(labelPsi,labelFis),function(x){ 
  return(c(mean(x[,"sd_psi"]),nrow(x)) ); 
}); 

colnames(eventCountResult)<-c("labelPsi","labelFis","sd_psi","count"); 
eventCountResult[,"count"]<-as.numeric(eventCountResult[,"count"]); 
#eventCountResult[] 

jpeg("psi_sd_fis.jpeg",width=900,height=600); 
p2<-ggplot(data=result2,aes(x=labelPsi,y=sd_psi,fill=factor(labelFis) ) ) + 
  geom_boxplot(position = position_dodge(width=0.9) )+theme_classic()+ 
geom_text(data = eventCountResult, aes(label = count),position = position_dodge(width=0.9) ); 

print(p2); 
dev.off(); 


lowFisResult<-result[result[,"labelFis"]=="fis<0.85",]; 
lowFisDensity<-ddply(lowFisResult,.(event_name),function(x){ 
  histVals<-hist(x[,"psi"],breaks=seq(0,1,0.01),xlim=c(0,1),plot=FALSE ); 
  density<-histVals$density; 
  
  return(c(x[1,"fis"],density,mean(x[,"psi"]) ) ); 
}); 

lowFisDensityOrder<-lowFisDensity[order(lowFisDensity[,102],decreasing = T),]; 


highFisResult<-result[result[,"labelFis"]=="fis>0.91",]; 
highFisDensity<-ddply(highFisResult,.(event_name),function(x){ 
  histVals<-hist(x[,"psi"],breaks=seq(0,1,0.01),xlim=c(0,1),plot=FALSE ); 
  density<-histVals$density; 
  
  return(c(x[1,"fis"],density,mean(x[,"psi"]) ) ); 
}); 

highFisDensityOrder<-highFisDensity[order(highFisDensity[,102],decreasing = T),]; 

pdf("fis_psi-density.pdf"); 
library(pheatmap); 

pheatmap(lowFisDensityOrder[,c(-1,-102)],cluster_rows=FALSE,
         color=colorRampPalette(c("aliceblue","darkred"))(50),
         cluster_cols=FALSE,show_rownames=F,show_colnames=T,
         main="low Fis score group(fis<0.85) order by PSI", 
         border_color=NA,labels_col=c(0,rep("",49),"","PSI",rep("",49),1) ); 

pheatmap(highFisDensityOrder[,c(-1,-102)],cluster_rows=FALSE,
         color=colorRampPalette(c("aliceblue","darkred"))(50),
         cluster_cols=FALSE,show_rownames=F,show_colnames=T,
         main="high Fis score group(fis>0.91) order by PSI",
         border_color=NA,labels_col=c(0,rep("",49),"","PSI",rep("",49),1) ); 
dev.off(); 

misoMergeRegion<-read.table("misoMergeRegion.csv",header=T,sep=",",as.is=T);
eventAnnova<-ddply(misoMergeRegion,.(event_name),function(x){
  
  pvalue<-NA; 
  tryCatch(pvalue<-oneway.test(psi ~V15 ,x)$p.value,error=function(e){return(NA);} ); 
  sample_count<-paste0(daply(x,.(V15),function(y){return(nrow(y));} ),collapse=",") 
  
  return(c(pvalue,sample_count)); 
});
eventAnnova[,2]<-as.numeric(eventAnnova[,2]);
write.table(eventAnnova,file="event_anova_p-value.tsv",sep="\t",row.names=F);

eventAnnovaNNA<-eventAnnova[!is.na(eventAnnova[,2]),]; 
eventAnnovaNNA[,2]<-as.numeric(eventAnnovaNNA[,2])
eventNamePValueList<-eventAnnovaNNA[order(eventAnnovaNNA[,2],decreasing = T)[1:20],1:2]; 

eventNameList<-eventNameList[,1]; 

miso_diseaseData<-read.table("miso_annodesp.predict.tsv",as.is=T,header=T); 

eventAnovaDiseaseData<-inner_join(eventAnnova,miso_diseaseData,by=c("event_name"="input")); 
write.table(eventAnovaDiseaseData[,c("event_name","V1","disease_probability")],
            file="event_anova_fis.tsv",sep="\t",row.names=F); 

eventAnovaDiseaseDataOrder<-eventAnovaDiseaseData[order(eventAnovaDiseaseData[,"V1"],decreasing = F),]; 

cor(eventAnovaDiseaseDataOrder[1:10,"V1"],eventAnovaDiseaseDataOrder[1:10,"disease_probability"]); 

pdf("event_box.pdf"); 
for(oneEventName in eventNameList){
  
misoMergeRegionEvent<-misoMergeRegion[misoMergeRegion[,"event_name"]==oneEventName,]; 

p<-ggplot(misoMergeRegionEvent)+
  geom_boxplot(aes(x=factor(V15),y=psi ) )+
  theme(axis.text.x = element_text(angle = 90, hjust = 1) )+ggtitle(substr(oneEventName,1,20 ))#+theme_classic(); 
print(p)

}
dev.off();
