
plotDPValue<-function(){
  
  phylopDataset<-alldatasetNA[,c("label","phylop")];
  ssDataset<-alldatasetNA[,c("label","ss_1","ss_2","ss_3","ss_4","ss_5","ss_6",
                           "ss_7","ss_8","ss_9","ss_10","ss_11","ss_12")];
  asaDataset<-alldatasetNA[,c("label","asa_1","asa_2","asa_3") ];
  disorderDataset<-alldatasetNA[,c("label","disorder_1","disorder_2",
                                 "disorder_3","disorder_4","disorder_5",
                                 "disorder_6","disorder_7","disorder_8",
                                 "disorder_9","disorder_10",
                                 "disorder_11","disorder_12")];

  pfamDataset<-alldatasetNA[,c("label","pfam1","pfam2")];
  ptmDataset<-alldatasetNA[,c("label","ptm")];
  
  features<-list();
  features[["phylop"]]<-phylopDataset;
  features[["ss"]]<-ssDataset;
  features[["asa"]]<-asaDataset;
  features[["disorder"]]<-disorderDataset;
  features[["pfam"]]<-pfamDataset;
  features[["ptm"]]<-ptmDataset;


  ksTestData<-data.frame(p=1:30, d=1:30,featureClass=1:30, featureName=1:30,p_raw=1:30 );
  
  featureIndex<-1;
  for(i in 1:length(features) ){
    currentFeatures<-features[[i]];
    pValues<-c();
  
    for(j in 2:ncol(currentFeatures)){
      p_raw<-ks.test(currentFeatures[currentFeatures[,1]=="NEUTRAL",j],
               currentFeatures[currentFeatures[,1]=="HGMD",j])$p.value;
      d<-ks.test(currentFeatures[currentFeatures[,1]=="NEUTRAL",j],
               currentFeatures[currentFeatures[,1]=="HGMD",j])$statistic;
      print(p);
      featureClass<-names(features)[i];
      featureName<-colnames(currentFeatures)[j];
    
      p<- -1*log10(p_raw);
      if(p>16)
        p<-16;
    
      ksTestData[featureIndex,]<-c( signif(p,3),signif(d,3),featureClass,featureName,p_raw );
      featureIndex<-featureIndex+1
    
      if(p!=0){
        p<- (-1)*log10(p)
      }else{
        p<-16
      }
      
      p<-rnorm(1,p,0.4);
      pValues<-c(pValues,p);
    }
    
  }

  ksTestData[,"p"]<-as.numeric(ksTestData[,"p"]); 
  ksTestData[,"d"]<-as.numeric(ksTestData[,"d"]); 

  colnames(ksTestData); 
  pdf("Each_feature_d_p_value.pdf"); 
  p1<-ggplot(ksTestData)+geom_point(aes(x=d,y=p,shape=featureClass))+
  #  geom_text(aes(x=d,y=p,label=featureName) )+
  theme_classic()+xlab("d-value")+ylab("p-value(-log10)")+
  ylim(c(0,17)) ; 
  print(p1);
  dev.off();
  write.table(ksTestData,file="ks-test.csv",sep=",",quote=FALSE,row.names=FALSE)
  
}

ksTestDataDesc<-cbind(ksTestData,abv(ksTestData[,"featureName"]) );

colnames(ksTestDataDesc)<-c("p-value(-log10)","D-value",
                            "feature class","feature symbol","p-value","feature description");
write.table(ksTestDataDesc,file="ks-test.tsv",sep="\t",quote=FALSE,row.names=FALSE);


plotasadisorder<-function(){
  
  pdf("asa_disorder.pdf");
  p1<-ggplot(alldatasetNA)+geom_point(aes(x=asa_1,y=disorder_2,color=label),shape=1)+theme_classic()+
  xlab("asa_1 (average ASA in translated amino acid sequence)" )+
  ylab("disorder_2 (max disorder score in amino acid sequence)" )+
  scale_colour_brewer(palette="Set1")+scale_color_identity();
  
  print(p1);
  dev.off();
  
}

plotcoildisorder<-function(){
  
  pdf("random_coil_disorder.pdf");
  p1<-ggplot(alldatasetNA)+geom_point(aes(x=ss_8,y=disorder_2,color=label),shape=1)+theme_classic()+
  xlab("ss_8 (min  probability of the amino acid in coil)")+
  ylab("disorder_2 (max disorder score in amino acid sequence)")+
  scale_colour_brewer(palette="Set1");
  
  print(p1);
  dev.off();
  
}

plotcoilasa<-function(){
  
  pdf("random_coil_asa.pdf");
  p1<-ggplot(alldatasetNA)+geom_point(aes(x=ss_8,y=asa_1,color=label),shape=1)+theme_classic()+
  xlab("ss_8 (min probability of the amino acid in coil)")+
  ylab("asa_1 (average ASA in translated amino acid sequence)")+
  scale_colour_brewer(palette="Set1")+scale_color_identity();
  
  print(p1);
  dev.off();
  
}

plotcoilphylop<-function(){
  
  pdf("random_coil_phylop.pdf");
  p1<-ggplot(alldatasetNA)+geom_point(aes(x=ss_8,y=phylop,color=label),alpha=0.3)+theme_classic()+
    xlab("ss_8 (min probability of the amino acid in coil)")+
    ylab("average phylop score")+
    scale_colour_brewer(palette="Set1");
  
  print(p1);
  dev.off();
  
}

plotcoilasa<-function(){
  
  pdf("random_coil_asa.pdf",width=10,height=10);
  coilDataset<-alldatasetNA[-1872,];
  rownames(coilDataset)<-paste0(sapply(str_split(rownames(coilDataset),"\\$"),"[",2),
                                sapply(str_split(rownames(coilDataset),"\\$"),"[",3),
                                sapply(str_split(rownames(coilDataset),"\\$"),"[",4))
  
  p1<-ggplot(coilDataset)+geom_point(aes(x=ss_8,y=asa_1,color=label),shape=1)+theme_classic()+
    xlab("ss_8 (min probability of the amino acid in coil)")+
    ylab("asa_1 (average ASA in translated amino acid sequence)")+
    geom_text(aes(x=ss_8,y=asa_1,label=rownames(coilDataset)),size=0.2 )+
    scale_colour_brewer(palette="Set1")#+scale_color_identity();
  
  print(p1);
  dev.off();
}
plotcoilasa();

plotcoilasaheat<-function(){
  pdf("random_coil_asa_heat.pdf");
  
  neutralData<-alldatasetNA[alldatasetNA[,"label"]=="NEUTRAL",];
  hgmdData<-alldatasetNA[alldatasetNA[,"label"]=="HGMD",];
  
  par(mfrow=c(1,2) );
  smoothScatter(neutralData[,"ss_8"],neutralData[,"asa_1"],
                xlab="ss_8 (min probability of the amino acid in coil)",
                  ylab="asa_1 (average ASA in translated amino acid sequence)",
  xlim=c(0,1),ylim=c(0,60) );
  
  smoothScatter(hgmdData[,"ss_8"],hgmdData[,"asa_1"],xlab="ss_8 (min probability of the amino acid in coil)",
                  ylab="asa_1 (average ASA in translated amino acid sequence)",
                  xlim=c(0,1),ylim=c(0,60) );
  
  dev.off();
}

plotbetasheetdisorder<-function(){
  pdf("beta-sheet_disorder.pdf");
  p1<-ggplot(alldatasetNA)+geom_point(aes(x=ss_5,y=disorder_2,color=label),shape=1)+theme_classic()+
    xlab("ss_5 (min probability of the amino acid in beta sheet structure)")+
    ylab("disorder_2 (max disorder score in amino acid sequence)")+
    scale_colour_brewer(palette="Set1");
  
  print(p1);
  dev.off();
  
  
}

plotalphahelixdisorder<-function(){
  pdf("alpha-helix_disorder.pdf");
  p1<-ggplot(alldatasetNA)+geom_point(aes(x=ss_10,y=disorder_2,color=label),shape=1)+theme_classic()+
    xlab("average probability of the amino acid in alpha-helix structure")+
    ylab("disorder_2 (max disorder score in amino acid sequence)")+
    scale_colour_brewer(palette="Set1");
  
  print(p1);
  #smoothScatter(alldatasetNA[,"ss_10"],alldatasetNA[,"disorder_2"]);
  dev.off();
  
}

plotbetasheetasa<-function(){
  
  jpeg("beta-sheet_alpha-helix_asa.jpeg",width=800,height=600);
  p1<-ggplot(alldatasetNA)+geom_point(aes(x=ss_4,y=asa_1,color=label),shape=1)+theme_classic()+
    xlab("ss_4 (average probability of the amino acid in beta sheet)")+
    ylab("asa_1 (average ASA in translated amino acid sequence)")+
    ggtitle("beta-sheet asa")+
    scale_colour_brewer(palette="Set1")#+scale_color_identity();
  
  p2<-ggplot(alldatasetNA)+geom_point(aes(x=ss_10,y=asa_1,color=label),shape=1)+theme_classic()+
    xlab("ss_10 (average probability of the amino acid in alpha helix)")+
    ylab("asa_1 (average ASA in translated amino acid sequence)")+
    ggtitle("alpha-helix asa")+
    scale_colour_brewer(palette="Set1")#+scale_color_identity();
  multiplot(p1,p2,cols=2);
    
  dev.off();
}
plotbetasheetasa();


library(rgl);

plotSecondaryStructure<-function(){
  pdf("Secondary structure figure.pdf");
  
  colors<-c();
  colors[which(alldatasetNA[,"label"]=="HGMD")]<-"red";
  colors[which(alldatasetNA[,"label"]=="NEUTRAL")]<-"blue";
  
  
  open3d();
  plot3d(alldatasetNA[,"ss_4"],alldatasetNA[,"ss_7"],alldatasetNA[,"ss_10"],
                xlab="average probability of beta sheet",
                ylab="average probability of random coil",
                zlab="average probability of alpha helix",
                col=colors);
  
  #dev.off();
}

plotcoilasa();
plotcoildisorder();
plotasadisorder();
plotcoilphylop();
plotcoilasaheat();
plotbetasheetdisorder();
plotalphahelixdisorder();

meASAHGMD<-median(alldataset[alldataset[,1]=="HGMD","asa_1"])
meASANeutral<-median(alldataset[alldataset[,1]=="NEUTRAL","asa_1"])

meASAHGMD<-35
meASANeutral<-35

sum(alldataset[alldataset[,1]=="HGMD" & (alldataset[,"asa_1"]>meASAHGMD),"ptm"]==0)
sum(alldataset[alldataset[,1]=="HGMD" & (alldataset[,"asa_1"]>meASAHGMD),"ptm"]>0)

sum(alldataset[alldataset[,1]=="HGMD" & (alldataset[,"asa_1"]<meASAHGMD),"ptm"]==0)
sum(alldataset[alldataset[,1]=="HGMD" & (alldataset[,"asa_1"]<meASAHGMD),"ptm"]>0)

sum(alldataset[alldataset[,1]=="NEUTRAL" & (alldataset[,"asa_1"]>meASANeutral),"ptm"]==0)
sum(alldataset[alldataset[,1]=="NEUTRAL" & (alldataset[,"asa_1"]>meASANeutral),"ptm"]>0)

sum(alldataset[alldataset[,1]=="NEUTRAL" & (alldataset[,"asa_1"]<meASANeutral),"ptm"]==0);
sum(alldataset[alldataset[,1]=="NEUTRAL" & (alldataset[,"asa_1"]<meASANeutral),"ptm"]>0);

scatterPtmHGMD<-alldataset[alldataset[,"label"]=="HGMD",c("asa_ave","ptm")];
scatterPtmNEUTRAL<-alldataset[alldataset[,"label"]=="NEUTRAL",c("asa_ave","ptm")];

#scatterPtmHGMD[scatterPtmHGMD[,"ptm"]>0,"ptm"]<-1;
#scatterPtmNEUTRAL[scatterPtmNEUTRAL[,"ptm"]>0,"ptm"]<-1;

#scatterPtmHGMD[,"ptm"]<- sapply( scatterPtmHGMD[,"ptm"],function(x){return(rnorm(1,x,0.000000001))}  );
#scatterPtmNEUTRAL[,"ptm"]<- sapply( scatterPtmNEUTRAL[,"ptm"],function(x){return(rnorm(1,x,0.000000001))} );

jpeg("asa_ptm.jpeg",width=1000,height=500);

a<-qplot(x=scatterPtmHGMD[,"ptm"],y=scatterPtmHGMD[,"asa_ave"],axes=T,
     ylab="average ASA",xlab="ptm",xlim=c(0,40),main="HGMD")#+scale_x_log10();

b<-qplot(x=scatterPtmNEUTRAL[,"ptm"],y=scatterPtmNEUTRAL[,"asa_ave"],axes=T,
     ylab="average ASA",xlab="ptm",xlim=c(0,40),main="NEUTRAL")#+scale_x_log10();

multiplot(a,b,cols=2);

dev.off();

jpeg("ptm_asa_scatter.jpeg");
d<-ggplot(alldataset)+geom_point( aes(x=ptm,y=asa_1,color=label));
print(d)
dev.off();


jpeg("ptm_asa_logic.jpeg");
alldatasetNAlogic<-alldatasetNA
#alldatasetNAlogic<-alldatasetNA[!is.na(alldatasetNA[,"ptm"]),];
alldatasetNAlogic[,"label"]<-factor(alldatasetNAlogic[,"label"],levels=c("NEUTRAL","HGMD") );
  
ptm_asa_logic<-glm(label~asa_1+ptm,alldatasetNAlogic,family="binomial");

#ptm_asa_logic<-glm(label~asa_1,alldatasetNAlogic,family="binomial");
ptm_asa_logic$residuals<-ptm_asa_logic$fitted.values-(as.numeric(ptm_asa_logic$data[,"label"])-1)

plot(ptm_asa_logic$data[names(ptm_asa_logic$residuals),"disorder_12"],ptm_asa_logic$residuals,
     ylab="residues",xlab="ptm_value" ,
     main=paste0(names(ptm_asa_logic$coefficients)[1],":",signif(ptm_asa_logic$coefficients[1],3),"\n",
                 names(ptm_asa_logic$coefficients)[2],":",signif(ptm_asa_logic$coefficients[2],3),"\n"
                 #names(ptm_asa_logic$coefficients)[3],":",signif(ptm_asa_logic$coefficients[3],3)
                 )
     )

dev.off();

asaMax<-max(alldataset[,"asa_ave"]);
asaMin<-min(alldataset[,"asa_ave"]);


ptmAsaFeature<-data.frame(hgmdPer=0,class="",label="",count=0,stringsAsFactors = F);

index<-1;
regions<-seq(from=asaMin,to=asaMax,by=0.5); 
#regions<-c(regions,62);
#regions<-c(0,50,62)
for(i in 2:length(regions)){
  hgmdRegion<-subset(alldataset,asa_ave>regions[i-1]&asa_ave<regions[i]&label=="HGMD");
  neutralRegion<-subset(alldataset,asa_ave>regions[i-1]&asa_ave<regions[i]&label=="NEUTRAL");
  
  hgmdPtmCount<-nrow(hgmdRegion);
  neutralPtmCount<-nrow(neutralRegion);
  
  hgmdPtmB0Count<-sum(hgmdRegion[,"ptm"]>0);
  neutralPtmB0Count<-sum(neutralRegion[,"ptm"]>0);
  
  ptmAsaFeature[index,]<-c(hgmdPtmB0Count/hgmdPtmCount,"HGMD",
                        paste0(signif(regions[i-1],2),"-",signif(regions[i],2) ),hgmdPtmCount);
  index<-index+1
  ptmAsaFeature[index,]<-c(neutralPtmB0Count/neutralPtmCount,"NEUTRAL",
                         paste0(signif(regions[i-1],2),"-",signif(regions[i],2)),neutralPtmCount );
  index<-index+1
  
}
options(digits=3) 

ptmAsaFeature[,"hgmdPer"]<-signif(as.numeric(ptmAsaFeature[,"hgmdPer"]),3)
ptmAsaFeature$label <- factor(ptmAsaFeature$label, levels = ptmAsaFeature$label);

p<-ggplot(ptmAsaFeature,aes(x=label, y= hgmdPer,fill=class ) ) + 
  geom_bar(position="dodge",stat = "identity")+
  ylab("PTM>0 percentile")+xlab("average asa in exon");

jpeg("asa_bin_ptm_percentile0-5.jpeg",width=1000,height=500);
print(p);
dev.off();
