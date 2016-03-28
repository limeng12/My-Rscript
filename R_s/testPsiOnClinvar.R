setwd(RPath);
source("getPSINeutralSNP.r");

predictClinvar5ByPsi<-function(){

  setwd(psiClinvar);
  snp<-read.table("clinvar_splicing_exon3_pathogenic.vcf",header=FALSE,as.is=TRUE);
  id<-snp[,3];
  
  setwd(psiClinvar);
  sortAndUniquePSI("clinvar_pathogenic_splicing_junction.tab");
  
  psiFileName<-"High dPSI unique.csv";
    
  psiFrame<-read.csv(file=psiFileName,header=TRUE,as.is=TRUE);
  colnames(psiFrame)<-c("id","transcript","dPSI","dPSI_percentile","PSI_WT");
  clnsig5Id<-intersect(id,psiFrame[,"id"]);
  
  psiFrame<-psiFrame[ psiFrame[,"id"] %in% clnsig5Id,];
  
  clinsig5Psi<-abs(psiFrame[,"dPSI"]);
  return(clinsig5Psi);
}

clinsig5PredictByPsi<-predictClinvar5ByPsi();

predictClinvar2ByPsi<-function(){
  
  setwd(psiClinvar);
  snp<-read.table("clinvar_splicing_exon3_benign.vcf",header=FALSE,as.is=TRUE);
  id<-snp[,3]
  
  setwd(psiClinvar);
  sortAndUniquePSI("clinvar_benign_splicing_junction.tab");
  
  psiFileName<-"High dPSI unique.csv";
  
  psiFrame<-read.csv(file=psiFileName,header=TRUE,as.is=TRUE);
  colnames(psiFrame)<-c("id","transcript","dPSI","dPSI_percentile","PSI_WT");
  clnsig2Id<-intersect(id,psiFrame[,"id"]);
  
  psiFrame<-psiFrame[ psiFrame[,"id"] %in% clnsig2Id,];
  
  clinsig2Psi<-abs(psiFrame[,"dPSI"]);
  return(clinsig2Psi);
  #hist(abs(psiFrame[,"dPSI"]) ,breaks=50);
}

clinsig2PredictByPsi<-predictClinvar2ByPsi();

setwd(workingDir);

#pdf("figure8.pdf");
jpeg("figure8.jpeg");

par(mfcol=c(2,1) );

hist(clinsig5PredictByPsi,breaks=20,
     main="SPANR on clinvar pathogenic splicing junction",
     xlab="|dPSI|"
     );

hist(clinsig2PredictByPsi,breaks=20,
     main="SPANR on clinvar benign splicing junction",
     xlab="|dPSI|"
     );

dev.off();

