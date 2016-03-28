library(randomForest);
setwd(RPath);
source("Global.r")
source("getPSINeutralSNP.r");
#splicing site in benign
#http://tools.genes.toronto.edu/results/0a790658-cf79-48ab-a3d9-7e68c824c90e


alldataset<-loadData(FALSE);
setwd(psiClinvar);

psiFileName<-sortAndUniquePSI("PSIpsiResultLinksclinvar_intron.forPSI.txt");
psiIntronData<-read.csv(psiFileName,header=TRUE,as.is=T);

psiFileName<-sortAndUniquePSI("psiResult_clinvar_splicing.txt"); 
psiSplicingData<-read.csv(psiFileName,header=TRUE,as.is=T);

psiData<-unique( rbind(psiIntronData,psiSplicingData) ); 

benignId<-read.table("clinvar_intron_benign.vcf",sep="\t",header=F,as.is=T); 
pathogenicId<-read.table("clinvar_annovar_multianno_splicing_pathogenic.vcf",sep="\t",header=F,as.is=T); 
splicingRsIds<-c(  benignId[,"V3"],pathogenicId[,"V3"] ); 
allSnpData<-rbind(benignId,pathogenicId); allSnpData[,1]<-paste0("chr",allSnpData[,1])

psiData<-psiData[is.element(psiData[,"id"],splicingRsIds),]

synonymous<-data.frame(
  id<-c("chr16:15815456:rs8046180:G-C","chr8:90970988:rs121908974:G-T","chr2:26455127:rs11552518:G-C"),
                       transcript<-c("chr16:MYH11:-:NM_022844:Exon32","chr8:NBN:-:NM_002485:Exon9","chr2:HADHA:-:NM_000182:Exon6"),
                       dPSI<-c(-32.59,-70.96,-73.31),
                       dPSI_percentile<-c(0.04,0.01,0.01),
                       PSI_WT<-c(81.05,81.8,70.24)
); 

colnames(synonymous)<-c("id","transcript","dPSI","dPSI_percentile","PSI_WT")
psiData<-rbind(psiData,synonymous)


transcriptId<-sapply(strsplit(psiData[,"transcript"] ,":"),"[",4); 
exonIndex<-substring(sapply(strsplit(psiData[,"transcript"] ,":"),"[",5),5 ); 
exonId<-paste0(transcriptId,"_",exonIndex);


psiDataAll<-cbind(psiData,exonId); 
psiDataAll[,"exonId"]<-as.character(psiDataAll[,"exonId"]); 

clinvarSplicing<-read.table("clinvar_splicing_psi",sep="\t",header=F,as.is=T); 
clinvarSplicing<-inner_join(clinvarSplicing,allSnpData,
                                by=c("V1"="V1","V2"="V2","V3"="V4","V4"="V5") ); 

exonId<-paste0(clinvarSplicing[,"V9"],"_",clinvarSplicing[,"V10"] ); 

clinvarSplicingData<-
  data.frame(id=".",transcript=clinvarSplicing[,"V9"],
             dPSI=clinvarSplicing[,5],dPSI_percentile=".",
             PSI_WT=".",exonId=exonId); 


psiDataAll<-rbind(psiDataAll,clinvarSplicingData); 

psiDataAll<-psiDataAll[order(abs(psiDataAll[,"dPSI"]),decreasing=T),]
psiDataAll<-psiDataAll[!duplicated(psiDataAll[,"exonId"]),]
cat(psiDataAll[,"exonId"],file=paste0(workingDir,"clinvarExonId"),sep="\n"); 

testUsingClinvar<-function(path,desp){
  
  setwd(path);
  #setwd(clinvarFeature2);
  
  despasass<-read.table(paste0(desp,".asa.ss.max_prob.data"),sep="\t",header=FALSE,as.is=TRUE);
  colnames(despasass)<-c("snpId","ss_1","ss_2","ss_3","ss_4","ss_5",
                         "ss_6","ss_7","ss_8","ss_9","ss_10","ss_11","ss_12",
                         "asa_1","asa_2","asa_3","exonId","ave_asa");
  
  despdisorder<-read.table(paste0(desp,".snp.ref.mut.seq.disorder-matrix"),sep="\t",header=FALSE,as.is=TRUE);
  colnames(despdisorder)<-c("snpId","region","disorder_1","disorder_2","disorder_3","disorder_4",
                            "disorder_5","disorder_6","disorder_7","disorder_8",
                            "disorder_9","disorder_10","disorder_11","disorder_12","exonId");
  
  desppfam<-read.table(paste0(desp,".snp.ref.mut.seq.pfam-matrix"),sep="\t",header=TRUE,as.is=TRUE);
  desppfam<-data.frame(id=desppfam[,"id"],pfam1=desppfam[,"overlap"],
                       pfam2=desppfam[,"coverage"],
                       exonId=desppfam[,"exonId"],stringsAsFactors = FALSE ); 
  
  despptm<-read.table(paste0(desp,".snp.ref.mut.seq.ptm-matrix.dbPTM"),sep="\t",header=TRUE,as.is=TRUE,quote=NULL);
  despptm<-data.frame(id=despptm[,"id"],ptm=rowSums(despptm[,3:(ncol(despptm)-1) ] ),
                      exonId=despptm[,"exonId"] ,stringsAsFactors = FALSE); 
  
  despPhylop<-read.table(paste0(desp,".snp.ref.mut.seq.mutBed.phylop"),sep="\t",header=FALSE,as.is=TRUE);
  colnames(despPhylop)<-c("snpId","phylop","r1","r2","r3","r4","exonId");
  
  proteinStarts<-as.numeric(sapply(strsplit(despdisorder[,"region"],":"),"[",2) );
  proteinEnds<-as.numeric(sapply(strsplit(despdisorder[,"region"],":"),"[",3) );
  proteinLength<-proteinEnds-proteinStarts+1
  despdisorder<-cbind(despdisorder,proteinLength);
  
  despData<-inner_join(despasass,despdisorder,by=c("snpId"="snpId","exonId"="exonId") );
  despData<-inner_join(despData,desppfam,by=c("snpId"="id","exonId"="exonId") );
  despData<-inner_join(despData,despptm,by=c("snpId"="id","exonId"="exonId") );
  despData<-inner_join(despData,despPhylop,by=c("snpId"="snpId","exonId"="exonId") );
  
  
  despData<-despData[!duplicated(despData[,"exonId"]),]; 
  rownames(despData)<-paste0(despData[,"snpId"],"$",despData[,"exonId"],"$",despData[,"region"],"$","desp"); 
  
  allFeatureNames<-c("phylop",
                     "ss_1","ss_2","ss_3","ss_4","ss_5",
                     "ss_6","ss_7","ss_8","ss_9","ss_10","ss_11","ss_12",
                     "asa_1","asa_2","asa_3",
                     "disorder_1","disorder_2","disorder_3","disorder_4",
                     "disorder_5","disorder_6","disorder_7","disorder_8",
                     "disorder_9","disorder_10","disorder_11","disorder_12",
                     "pfam1","pfam2","ptm","proteinLength"); 
  
  
  if(path==clinvarFeature2PSI20){
    
    outTable<-inner_join(despData,psiDataAll,by=c("exonId"="exonId") ); 
    exonIdTrain<-readLines(paste0(workingDir,"training_exonId.txt")); 
    inTraining<-sapply( outTable[,"exonId"],is.element,exonIdTrain ); 
    outTable<-cbind(outTable,inTraining);
  
    despData<-despData[!is.element(despData[,"exonId"],exonIdTrain),]; 
    despData<-outTable[!outTable[,"inTraining"],]
    despData<-outTable;
    despData<-despData[,allFeatureNames]; 
  
    predictClinvar<-predict(modelrandomforestAll,despData,type="prob")[,1]; 
    outTable<-cbind(outTable,predictClinvar); 
  
    outTable<-outTable[abs(outTable[,"dPSI"])>10,]; 
    write.table(outTable[,c("exonId","predictClinvar","dPSI","inTraining",allFeatureNames) ],
              file="result_table.tsv",sep="\t",row.names = F,quote=F); 
    #return(outTable[!outTable[,"inTraining"],"predictClinvar"]);
   return(outTable[,"predictClinvar"])
  }else{
    outTable<-inner_join(despData,psiDataAll,by=c("exonId"="exonId") ); 
    outTable<-outTable[abs(outTable[,"dPSI"])>10,]; 
    despData<-outTable[,allFeatureNames];
    
    predictClinvar<-predict(modelrandomforestAll,despData,type="prob")[,1]; 
    
    data<-cbind(predictClinvar,despData);
    write.table(data,file="pahtogenic.tsv",sep="\t");
    return(predictClinvar);
  }
  #return(predictClinvar);
}

setwd(workingDir);
#pdf("figure7.pdf");

#predictClinvar2<-testUsingClinvar(clinvarFeature2);
#predictClinvar5<-testUsingClinvar(clinvarFeature5,"clinvardesp");
predictClinvar5SplicingJunction<-testUsingClinvar(clinvarFeature5SplicingJunction,"clinvardesp");
predictClinvar2SplicingJunction<-testUsingClinvar(clinvarFeature2PSI0,"clinvar_benign_psi0desp");

predictClinvar2<-predictClinvar2SplicingJunction;
predictClinvar5<-predictClinvar5SplicingJunction;

setwd(workingDir);

pdf("figure_5_clinvar_test.pdf");
p1<-ggplot( )+geom_histogram(aes(x=predictClinvar2,y=..density..),binwidth =0.05)+
  ggtitle("pathogenic in splicing junction")+
  #geom_density(aes(x=predictClinvar2),adjust=0.4)+
  xlab("disease causing probability")+ggtitle("Predicted on pathogenic splicing junction");

p2<-ggplot( )+geom_histogram(aes(x=predictClinvar5,y=..density..),binwidth =0.05)+      
  #geom_density(aes(x=predictClinvar5),adjust=0.4)+
  ggtitle("Predicted on Clinvar pathogenic")+xlab("disease causing probability");


#print(p2);
multiplot(p1,p2,cols=1 );
dev.off();

pdf("figure6_test_on_clinvar.pdf");

ksPvalue<-ks.test(predictClinvar2,predictClinvar5)$p.value ;
par(mfrow=c(2,1));

hist(predictClinvar2,main="Predicted on Clinvar benign ( |dPSI|>20%)", 
     breaks=10,xlab="disease causing probability",col="gray");

hist(predictClinvar5,main="Predicted on Clinvar pathogenic (splicing junction)" ,
     breaks=10,xlab="disease causing probability",col="gray");

#print(p2);
#multiplot(p1,p2,cols=1 );


labelBenign<-rep("NEUTRAL",length(predictClinvar2) );
labelPathogenic<-rep("HGMD",length(predictClinvar5) );

rocData<-data.frame(value=c(predictClinvar2,predictClinvar5),label=c(labelBenign,labelPathogenic) );

tprFun<-function(x){
  t<-sum(rocData[rocData$label=="HGMD","value"]>x)/sum(rocData$label=="HGMD");
  return(t);
}

fprFun<-function(x){
  t<-sum(rocData[rocData$label=="NEUTRAL","value"]>x)/sum(rocData$label=="NEUTRAL");
  return(t);
}

thresholds<-seq(-0.1,1.1,0.001);
tpr<-sapply(thresholds,tprFun);
fpr<-sapply(thresholds,fprFun);

aucValue<-mean(sample(rocData[rocData$label=="HGMD","value"],1000000,replace=T) >
                 sample(rocData[rocData$label=="NEUTRAL","value"],1000000,replace=T) );
dev.off();
sum(predictClinvar5>0.9)/length(predictClinvar5);

