#Global variables

RPath<-"/Users/mengli/Documents/projects/splicingSNP/R/";

HGMDFeaturePath<-"/Users/mengli/Documents/projects/splicingSNP/features/hgmd";
HGMDEnsemblFeaturePath<-"/Users/mengli/Documents/projects/splicingSNP/features/ensembl_hgmd/";

yaoqiNeutralFeaturePath<-"/Users/mengli/Documents/projects/splicingSNP/features/normalYaoqi/";
yaoqiNeutralEnsemblFeaturePath<-"/Users/mengli/Documents/projects/splicingSNP/features/ensembl_normal/";

yaoNeutralPath<-"/Users/mengli/Documents/projects/splicingSNP/features/normalYaoqi/";
psiNeutralFeaturePath<-"/Users/mengli/Documents/projects/splicingSNP/features/1000genomePSI/AF5population/";

psiFilesDir<-"/Users/mengli/Documents/projects/splicingSNP/features/1000genomePSI/AF5population/";

psiFileName<-"sorted psi file name";
workingDir<-"/Users/mengli/Documents/projects/splicingSNP/features/YaoqiExonshgmd/";


refGeneBedFilePath<-"/Users/mengli/Documents/projects/splicingSNP/annotation/refgene.bed";

hgmdRawDataPath<-"/Users/mengli/Documents/projects/splicingSNP/hgmdRawData";
hgmdSnvSpliceJunc<-"/Users/mengli/Documents/projects/splicingSNP/hgmdRawData/hgmdSplicingJunction/";
hgmdSnvSpliceJunc25<-"/Users/mengli/Documents/projects/splicingSNP/hgmdRawData/hgmdSplicingJunction25/";


#clinvarFeature<-"/Users/mengli/Documents/projects/splicingSNP/code/clinvarFeatures/"; 
#clinvarFeature1<-"/Users/mengli/Documents/projects/splicingSNP/code/clinvarFeatures1/"; 
#clinvarFeature2<-"/Users/mengli/Documents/projects/splicingSNP/features/newpsiClinvar/clinvar_benign_features_psi5"; 
#clinvarFeature5<-"/Users/mengli/Documents/projects/splicingSNP/features/newpsiClinvar/clinvar_pathogenic_features_psi5/"; 
clinvarFeature5SplicingJunction<-"/Users/mengli/Documents/projects/splicingSNP/features/clinvar5SJ/"; 
clinvarFeature2PSI20<-"/Users/mengli/Documents/projects/splicingSNP/features/clinvar2psi20"; 
clinvarFeature2PSI0<-"/Users/mengli/Documents/projects/splicingSNP/features/clinvar2psi0"; 

psiClinvar<-"/Users/mengli/Documents/projects/splicingSNP/features/psiClinvar/"; 

clinVarPath<-paste(psiClinvar,"clinvar_splicing_exon3.vcf",sep=""); 

normal1000GenomeSplicingPSI5Features<-"/Users/mengli/Documents/projects/splicingSNP/code/normal1000GenomeSplicingPSI5Features/"; 
normal1000GenomeSplicingPSI20Features<-"/Users/mengli/Documents/projects/splicingSNP/code/normal1000GenomeSplicingPSI20Features/"; 


splicingJunction1000Genome<-"/Users/mengli/Documents/projects/splicingSNP/1000genome/SplicingJunction2/";

#all the data set, 6 classes, 30 features.
alldataset<-NULL;

#all the data set, 6 classes, 30 features, missing ptm values are replaced by NA.
alldatasetNA<-NULL;

#2/3 of all the data set.
allDatasetTrain<-NULL;

#1/3 of all the data set.
allDatasetTest<-NULL;

#feature selection process and their AUC.
featureAUC<-NULL;

#the distance of HGMD snp's distance to a exon.
threshold<-NULL;#HGMD max distance to a exon.

#predict on clinvar data set.
predictClinvar<-NULL;

#predict on PSI data set.
predictPSI<-NULL;

clinsig5PredictByPsi<-NULL;
clinsig2PredictByPsi<-NULL;
importantFeatures<-NULL;

setwd("/Users/mengli/Documents/projects/splicingSNP/R");
#abbreviation

#hgmd gene ids
hgmdIds<-NULL;

source("abv.r")


