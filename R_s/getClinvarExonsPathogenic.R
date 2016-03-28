library('BSgenome.Hsapiens.UCSC.hg19')
setwd(RPath);
source("getPSINeutralSNP.r")
source("Global.r");
source("BED.r")
source("VCF.R");

#Total 111000
#splicing junction 2228
#clnSig=1 617
#clnSig=2 75
#clnSig=3 6
#clnSig=4 331
#clnsig=5 1108
#clnsig=255 48

getClinvarExons<-function(clinvarFile){
  setwd(psiClinvar);
  
  vcf<-new("VCF");
  vcf<-readFromVCFFile(vcf,clinvarFile );
  #only keep the splice junction snp
  #vcf<-vcfFilterByAttr(vcf,"splicing");
  
  #benign and likely benign filter
  #vcf<-vcfFilterByAttr(vcf,"CLNSIG=3;|CLNSIG=2;");
  #disease filter
  #vcf<-rmIndel(vcf);
  pattern5<-"CLNSIG=[^;]*[^25]{1}5[^;]*;|CLNSIG=5[^;]*;";
  pattern2<-"CLNSIG=[^;]*2[^5]{1}[^;]*;|CLNSIG=[^;]*2;";
  ##pattern test
  #2|2
  #2,2
  #2;
  #
  #CLNSIG=
  #3
  #5|2
  #5,2
  #
  
  vcf<-vcfFilterByAttr(vcf,pattern5);
  #"CLNSIG=[^;]*[^25]{1}5[^;]*;|CLNSIG=5;"
  #"CLNSIG=[^;]*2[^5]{1}[^;]*;"
  
  #"CLNSIG=5;|CLNSIG=[^;]*\\|5;|CLNSIG=[^;]*\\|5\\||CLNSIG=5\\|"
  #"CLNSIG=2;|CLNSIG=[^;]*\\|2;|CLNSIG=[^;]*\\|2\\||CLNSIG=2\\|"
  
  #remove indels
  #vcf<-rmIndel(vcf);
  vcfFrame<-getDataFrame(vcf);
  write.csv(vcfFrame,file="clinvar5_splicing_junction.vcf",quote=TRUE,row.names=FALSE);
  #CLNSIG
  
}

getClinvarExons(clinVarPath); 

##Intron SNP predict by SPANR 
setwd(psiClinvar);
psidata<-read.table("clinvar_splicing_psi",header=F,sep="\t",as.is=T); 

clinvarPathogenicData<-read.table("clinvar5_splicing_junction.vcf",header=TRUE,as.is=TRUE,sep=","); 
clinvarPathogenicData[,"chr"]<-paste0("chr",clinvarPathogenicData[,"chr"]); 

psidatapathogenic<-inner_join(clinvarPathogenicData,psidata,by=c("chr"="V1","pos"="V2","ref"="V3","alt"="V4") ); 

psidatapathogenicnFilename<-"sort_unique_psi_clinvar_pathogenic"; 
write.csv(psidatapathogenic,file=psidatapathogenicnFilename,row.names=FALSE); 

#psi20<-psidatapathogenic[abs(psidatapathogenic[,"V5"])>20,]; 
#write.csv(psi20,file="sort_unique_psi20_clinvar_pathogenic",row.names=FALSE); 

setwd(psiClinvar); 
getNeutralExonCentralPos2(psidatapathogenicnFilename,0,"clinvar_pathogenic_psi0"); 

