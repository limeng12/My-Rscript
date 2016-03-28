library('BSgenome.Hsapiens.UCSC.hg19')
setwd(RPath);
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
  vcf<-rmIndel(vcf);
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
  
  vcf<-vcfFilterByAttr(vcf,pattern2);
  #"CLNSIG=[^;]*[^25]{1}5[^;]*;|CLNSIG=5;"
  #"CLNSIG=[^;]*2[^5]{1}[^;]*;"
  
  #"CLNSIG=5;|CLNSIG=[^;]*\\|5;|CLNSIG=[^;]*\\|5\\||CLNSIG=5\\|"
  #"CLNSIG=2;|CLNSIG=[^;]*\\|2;|CLNSIG=[^;]*\\|2\\||CLNSIG=2\\|"
  
  #remove indels
  #vcf<-rmIndel(vcf);
  vcfFrame<-getDataFrame(vcf);
  write.csv(vcfFrame,file="filter.vcf",quote=FALSE,row.names=FALSE);
  #CLNSIG
  
  snp<-data.frame(chr=vcfFrame[,1],pos=vcfFrame[,2],strand=rep("*",nrow(vcfFrame) )  );
  write.csv(snp,file="snp.csv",row.names=FALSE,quote=FALSE);
  b=new("BED");
  mBed<-readFromBedFile(b,refGeneBedFilePath);
  
  #threshold here is the
  threshold<-10;
  
  clinvarExons<-getNearestExons(mBed,snp,threshold);
  clinvarExons<-unique(clinvarExons);
  
  chr<-clinvarExons[,"chrosome"];
  exonBegPos<-as.numeric(clinvarExons[,"exonBegPos"] );
  exonEndPos<-as.numeric(clinvarExons[,"exonEndPos"] );
  strand<-clinvarExons[,"strand"];
  
  exonCenterPos<-round( (exonBegPos+exonEndPos)/2  );
  refSeq<-as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19,chr,exonCenterPos,exonCenterPos) );
  
  exonRegions<-str_c(chr,":",exonBegPos,"-",exonEndPos);
  xinJunInput<-str_c(chr,":",exonCenterPos,"\t",refSeq,"\t","T","\t",strand);
  
  cat(exonRegions,file="clinvarOld",sep="\n");# for debuging
  cat(xinJunInput,file="clinvar2_splicing_junction",sep="\n");

}

getClinvarExons(clinVarPath); 
#getClinvarExons("/home/limeng/splicingSNP/features/psiClinvar/clinvar_annovar.hg19_multianno_splicing.vcf"); 


##Intron SNP predict by SPANR
setwd(psiClinvar);
psiFileName<-sortAndUniquePSI("PSIpsiResultLinksclinvar_intron.forPSI.txt"); 
#psiFileName="High dPSI unique.csv";
psidata<-read.csv(psiFileName,header=TRUE); 
rownames(psidata)<-psidata[,1]; 

clinvarBenignData<-read.table("clinvar_intron_benign.vcf",header=TRUE,as.is=TRUE); 
ids<-clinvarBenignData[,3];
overlapIds<-intersect(ids,psidata[,1]);
psidatabenign<-psidata[overlapIds,];
#
#chr16:15815456:rs8046180:G-C  chr16:MYH11:-:NM_022844:Exon32	-32.59	0.04	81.05
#chr8:90970988:rs121908974:G-T	chr8:NBN:-:NM_002485:Exon9	-70.96	0.01	81.8
#chr2:26455127:rs11552518:G-C	chr2:HADHA:-:NM_000182:Exon6	-73.31	0.01	70.24

#
#
synonymous<-data.frame(id<-c("chr16:15815456:rs8046180:G-C","chr8:90970988:rs121908974:G-T","chr2:26455127:rs11552518:G-C"),
                       transcript<-c("chr16:MYH11:-:NM_022844:Exon32","chr8:NBN:-:NM_002485:Exon9","chr2:HADHA:-:NM_000182:Exon6"),
                       dPSI<-c(-32.59,-70.96,-73.31),
                       dPSI_percentile<-c(0.04,0.01,0.01),
                       PSI_WT<-c(81.05,81.8,70.24)
                       ); 

colnames(synonymous)<-c("id","transcript","dPSI","dPSI_percentile","PSI_WT"); 

psidatabenign<-rbind(psidatabenign,synonymous); 

psidatabenignFilename<-"sort_unique_psi_clinvar_benign"; 
write.csv(psidatabenign,file=psidatabenignFilename,row.names=FALSE); 

psi20<-psidatabenign[abs(psidatabenign[,3])>0,]; 
write.csv(psi20,file="sort_unique_psi20_clinvar_benign",row.names=FALSE); 

setwd(psiClinvar); 
options(warn=2)
getNeutralExonCentralPos2(psidatabenignFilename,0,"clinvar_benign_psi20"); 
