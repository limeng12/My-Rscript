setwd(RPath);
source("getPSINeutralSNP.r");

setwd("/home/limeng/Projects/xinjun/");
data<-read.csv("benign.csv");

#20 is the psi threshold
getNeutralExonCentralPos2("benign.csv",10);

setwd("/home/limeng/splicingSNP/code/benignExonbodyPSI20/");

normalasass<-read.table("normal.asa.ss.max_prob.data",sep="\t",header=FALSE);
normalasass<-normalasass[!duplicated(normalasass[,1]),];rownames(normalasass)<-(normalasass[,1]);
normaldisorder<-read.table("normal.snp.ref.mut.seq.disorder-matrix",sep="\t",header=FALSE);
normaldisorder<-normaldisorder[!duplicated(normaldisorder[,1]),];rownames(normaldisorder)<-(normaldisorder[,1]);            
normalpfam<-read.table("normal.snp.ref.mut.seq.pfam-matrix",sep="\t",header=TRUE);
normalpfam<-normalpfam[!duplicated(normalpfam[,1]),];rownames(normalpfam)<-(normalpfam[,1]);
normalptm<-read.table("normal.snp.ref.mut.seq.ptm-matrix.dbPTM",sep="\t",header=TRUE,quote=NULL);
normalptm<-normalptm[!duplicated(normalptm[,1]),];rownames(normalptm)<-(normalptm[,1]);

normalPhylop<-read.table("normal.snp.ref.mut.seq.mutBed.phylop",sep="\t",header=FALSE);
normalPhylop<-normalPhylop[!duplicated(normalPhylop[,1]),];rownames(normalPhylop)<-(normalPhylop[,1]);


cat(setdiff(rownames(normalptm),rownames(normaldisorder)) ,file="wrong exons.txt",sep="\n" );
setwd(workingDir);


normalSNP<-intersect(normalasass[,1],normaldisorder[,1]);
normalSNP<-intersect(normalSNP,normalpfam[,1]);
normalSNP<-intersect(normalSNP,normalptm[,1]);
normalSNP<-intersect(normalSNP,normalPhylop[,1]);


normalData<-cbind(normalPhylop[normalSNP,2],normalasass[normalSNP,2:16],
                  normaldisorder[normalSNP,3:14],normalpfam[normalSNP,ncol(normalpfam)],
                  rowSums(normalptm[normalSNP,3:ncol(normalptm)]) );

colnames(normalData)<-c("phylop",
                        "ss_1","ss_2","ss_3","ss_4","ss_5",
                        "ss_6","ss_7","ss_8","ss_9","ss_10","ss_11","ss_12",
                        "asa_1","asa_2","asa_3",
                        "disorder_1","disorder_2","disorder_3","disorder_4",
                        "disorder_5","disorder_6","disorder_7","disorder_8",
                        "disorder_9","disorder_10","disorder_11","disorder_12",
                        "pfam", "ptm");

modelrandomforest<-randomForest(formula=label~., data=alldataset,
                                ntree=500,proximity=TRUE,
                                replace=FALSE,nodesize=20,maxnodes=70);

predictClinvarBenign<-predict(modelrandomforest,normalData,type="prob")[,1];


