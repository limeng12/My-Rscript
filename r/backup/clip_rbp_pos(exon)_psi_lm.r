library(dplyr)
library(stringr)
library(Matrix)
library(qlcMatrix)
library(pheatmap)

setwd("E:\\limeng\\RBP_network\\r");

clip_binding_data<-read.table("../anno/gene_clip_target_se_exon.tsv",sep="\t",as.is=T,header=T,row.names=1);
clip_binding_data[,"exonId"]<-str_sub(clip_binding_data[,"exonId"],6,-1);

clip_names<-colnames(clip_binding_data)[6:ncol(clip_binding_data)];

expression_psi_data<-read.table(file="../anno/all_data_expression_psi.tsv",as.is=T,header=T);
sample_names<-colnames(expression_psi_data);

expression_psi_data<-cbind(rownames(expression_psi_data),expression_psi_data);

colnames(expression_psi_data)[1]<-"name";
expression_psi_data[,"name"]<-as.character(expression_psi_data[,"name"]);

clip_binding_exp_psi_data<-inner_join(clip_binding_data,expression_psi_data,by=c("exonId"="name") );

write.table(clip_binding_exp_psi_data,file="../anno/clip_binding_psi.tsv",sep="\t",row.names=F,col.names=T);

median_psi<-apply(clip_binding_exp_psi_data[,sample_names ],1,function(x){
  return(median(as.numeric(x),na.rm=T) );
});

clip_binding_exp_psi_data<-cbind(clip_binding_exp_psi_data,median_psi);

clip_binding_median_psi<-clip_binding_exp_psi_data[,c(clip_names, "median_psi" )];
clip_binding_median_psi_lm<-lm(median_psi~.,data=clip_binding_median_psi);
clip_binding_median_psi_lm_sum<-summary(clip_binding_median_psi_lm);

clip_binding_median_psi_lm_pValue<--1*log10(clip_binding_median_psi_lm_sum$coefficients[,4]);

clip_binding_median_psi_lm_pValue<-sort(clip_binding_median_psi_lm_pValue)[1:(length(clip_binding_median_psi_lm_pValue)-1)];

pdf("../anno/linear regression of RBP(exon) and Event.pdf",width=8,height=5);
barplot(clip_binding_median_psi_lm_pValue,las=2,main="Linear regression P-value of RBP and Event(-log10)");

dev.off();

clip_binding<-clip_binding_median_psi[,-ncol(clip_binding_median_psi)];
#clip_binding_cor_mat<-t(as.matrix(clip_binding));
clip_binding_cor_mat<-Matrix(as.matrix(clip_binding) );

clip_binding_cor<-corSparse(clip_binding_cor_mat);

colnames(clip_binding_cor)<-colnames(clip_binding)
rownames(clip_binding_cor)<-colnames(clip_binding)

pdf("../anno/RBP(exon)_cor.pdf");
pheatmap(clip_binding_cor,cluster_rows=T,cluster_cols=T);
dev.off();

#clip_binding_cor_binary<-sapply(clip_binding_cor,function(x){if(x<0)return(-1);return(1); } );
#pdf("../anno/RBP(exon)_cor_binary.pdf");
#pheatmap(clip_binding_cor_binary,cluster_rows=F,cluster_cols=F);
#dev.off();
