library(dplyr)
library(stringr);

setwd("E:\\limeng\\RBP_network\\r");

clip_binding_data<-read.table("../anno/gene_clip_target_se.tsv",sep="\t",as.is=T,header=T,row.names=1);
clip_binding_data[,"exonId"]<-str_sub(clip_binding_data[,"exonId"],6,-1);

expression_psi_data<-read.table(file="../anno/all_data_expression_psi.tsv",as.is=T,header=T);

expression_psi_data<-cbind(rownames(expression_psi_data),expression_psi_data);

colnames(expression_psi_data)[1]<-"name";
expression_psi_data[,"name"]<-as.character(expression_psi_data[,"name"]);

clip_binding_exp_psi_data<-inner_join(clip_binding_data,expression_psi_data,by=c("exonId"="name") );

write.table(clip_binding_exp_psi_data,file="../anno/clip_binding_psi.tsv",sep="\t",row.names=F,col.names=T);

median_psi<-apply(clip_binding_exp_psi_data[,55:(ncol(clip_binding_exp_psi_data)-1) ],1,function(x){
  return(median(as.numeric(x),na.rm=T) );
});

clip_binding_exp_psi_data<-cbind(clip_binding_exp_psi_data,median_psi);

clip_binding_median_psi<-clip_binding_exp_psi_data[,c(6:53,169)];
clip_binding_median_psi_lm<-lm(median_psi~.,data=clip_binding_median_psi);
clip_binding_median_psi_lm_sum<-summary(clip_binding_median_psi_lm);

clip_binding_median_psi_lm_pValue<--1*log10(clip_binding_median_psi_lm_sum$coefficients[,4]);

clip_binding_median_psi_lm_pValue<-sort(clip_binding_median_psi_lm_pValue)[1:(length(clip_binding_median_psi_lm_pValue)-1)];

pdf("../anno/linear regression of RBP and Event.pdf",width=8,height=5);
barplot(clip_binding_median_psi_lm_pValue,las=2,main="Linear regression P-value of RBP and Event(-log10)");

dev.off();


clip_binding<-clip_binding_median_psi[,-ncol(clip_binding_median_psi)];
clip_binding_cor_mat<-t(as.matrix(clip_binding));
clip_binding_cor<-cor(clip_binding_cor_mat);


pdf("../anno/RBP_cor.pdf");
pheatmap(clip_binding_cor,cluster_rows=F,cluster_cols=F);
dev.off();

clip_binding_cor_binary<-sapply(clip_binding_cor,function(x){if(x<0)return(-1);return(1); } );
pdf("../anno/RBP_cor_binary.pdf");
pheatmap(clip_binding_cor_binary,cluster_rows=F,cluster_cols=F);
dev.off();
