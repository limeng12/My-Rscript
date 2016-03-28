#library("devtools");
#install_github("kassambara/factoextra");
library("factoextra");
#source("BED.r")

#install_github("limeng12/pca3d")
#library(pca3d)
setwd(RPath)
source("object3d.r")
source("pca3d.r")
source("helper.r")

plotFigure4<-function(taxes){
  
  #jpeg("figure3.jpeg",width=1000,height=800);
  
  par(cex.lab=0.5, cex.axis=0.5, cex.main=0.5, cex.sub=0.5);
  pcaAllDataSet<-alldataset;
  
  selectNames<-c("Label","phylop","ss_1","ss_2","ss_3",
                             "ss_4","ss_5","ss_6","ss_7","ss_8",
                             "ss_9","ss_10","ss_11","ss_12",
                             "asa_1","asa_2","asa_3",
                             "d_1","d_2","d_3","d_4","d_5","d_6","d_7",
                             "d_8","d_9","d_10","d_11","d_12",
                             "pfam1","pfam2",
                             "ptm","proteinLength"
  );
  #selectNames
  
  #rownames(pcaAllDataSet)<-NULL;
  #names<-colnames(pcaAllDataSet);
  
  a<-princomp(pcaAllDataSet[,-1],cor=TRUE);
  
  p<-fviz_pca_biplot(a,axes=taxes,data=pcaAllDataSet[,-1],label ="var",habillage=pcaAllDataSet[,1],
                  col.var = c(rep("blue",1),rep("red",12),
                              rep("darkgoldenrod2",3),rep("black",12),
                              rep("green",2),rep("blueviolet",1)  ,rep("yellow",1))
                  ,labelsize=2)+ theme(text = element_text(size=15) );
  p<-p+ggtitle("PCA biplot");
  
  
  p<-p+geom_point(aes(x=seq(0.0001,0.0001,by=0.0001),y=seq(0.0001,0.0001,by=0.0001),
                      size=c("phylop","secondary structure","ASA","disorder","pfam","ptm","protein Length") ),
                  colour="black")
  p<-p+guides(size=guide_legend("Features", 
                                override.aes=list(shape=95, size = 10,
                                colour=c("darkgoldenrod2","black","green","blue","yellow","blueviolet","red") )));
  print(p);
  return(p);
  
}

setwd(workingDir);

pdf("figure4_pca_correlation_analysis.pdf",width=12,height=10);
p12<-plotFigure4(c(1,2) );
#p13<-plotFigure4(c(1,3) );
#p14<-plotFigure4(c(1,4) );
#p23<-plotFigure4(c(2,3) );
#p24<-plotFigure4(c(2,4) );
#p34<-plotFigure4(c(3,4) );

#multiplot(p12,p13,p14,p23,p24,p34 ,cols=2);
dev.off();


plotFigure43d<-function(taxes){
  
  #jpeg("figure3.jpeg",width=1000,height=800);
  
  par(cex.lab=0.5, cex.axis=0.5, cex.main=0.5, cex.sub=0.5);
  pcaAllDataSet<-alldataset;
  
  colnames(pcaAllDataSet)<-c("label","phylop","ss_1","ss_2","ss_3",
                             "ss_4","ss_5","ss_6","ss_7","ss_8",
                             "ss_9","ss_10","ss_11","ss_12",
                             "asa_1","asa_2","asa_3",
                             "d_1","d_2","d_3","d_4","d_5","d_6","d_7",
                             "d_8","d_9","d_10","d_11","d_12",
                             "pfam1","pfam2",
                             "ptm"
  );
  
  #names<-colnames(pcaAllDataSet);
  palColors<-c("darkgoldenrod1","dodgerblue");
  colors<-palColors[factor(pcaAllDataSet[,"label"])];
  a<-princomp(pcaAllDataSet[,-1],cor=TRUE);
  pca3d(a,col=colors,biplot=TRUE,radius=0.25,biplot.vars=1:31,
        cex=0.5,loading_col=c(rep("blue",1),rep("red",12),
                              rep("darkgoldenrod2",3),rep("black",12),
                              rep("green",2),rep("blueviolet",1)  ) );
  
  snapshotPCA3d("pca-3d.png");
  rgl.postscript("pca-3d.pdf",fmt="pdf");
  
}


getCluster<-function(){
  pcaAllDataSet<-alldataset;
  colnames(pcaAllDataSet)<-c("Label","phylop","ss_1","ss_2","ss_3",
                             "ss_4","ss_5","ss_6","ss_7","ss_8",
                             "ss_9","ss_10","ss_11","ss_12",
                             "asa_1","asa_2","asa_3",
                             "d_1","d_2","d_3","d_4","d_5","d_6","d_7",
                             "d_8","d_9","d_10","d_11","d_12",
                             "pfam",
                             "ptm"
  );
  
  names<-colnames(pcaAllDataSet);
  a<-princomp(pcaAllDataSet[,-1],cor=TRUE);
  scores<-a$scores;
    
  smallCluster<-rownames(scores)[ (scores[,1]>4)&(scores[,1]<6)&(scores[,2]>-2)&(scores[,2]< (-1) )];
  smallCluster<-smallCluster[grep("HGMD",smallCluster)]
  
  snpStrs<-strsplit(smallCluster,"\\$");
  
  exonIds<-str_c(sapply(snpStrs,"[",2)," ",sapply(snpStrs,"[",3) );
  cat(exonIds,file="cluster_exon_ids.txt",sep="\n");
  
  transcriptId<-sapply(strsplit(sapply(snpStrs,"[",3),":") ,"[",1); 
  
  cat(unique(transcriptId),file="target_id.txt",sep="\n");
  #mBed=new("BED");
  #mBed<-readFromBedFile(mBed,"../../annotation/refgene.bed");
  #transcriptId_exon<-getTranscriptIdExonNumber(mBed,snps);
  
}

getCluster();
