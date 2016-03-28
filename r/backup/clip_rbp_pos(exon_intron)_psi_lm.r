library(dplyr)
library(stringr)
library(Matrix)
library(qlcMatrix)
library(pheatmap)
library(ggplot2)

setwd("/Users/mengli/Documents/projects/RBP_network/r");

clip_binding_data_exon<-read.table("../anno/gene_clip_target_se_exon.tsv",sep="\t",as.is=T,header=T,row.names=1);
clip_binding_data_exon[,"exonId"]<-str_sub(clip_binding_data_exon[,"exonId"],6,-1);

clip_names_exon<-paste0(colnames(clip_binding_data_exon)[6:ncol(clip_binding_data_exon)],"_exon");
colnames(clip_binding_data_exon)[6:ncol(clip_binding_data_exon)]<-clip_names_exon


clip_binding_data_intron_up<-read.table("../anno/gene_clip_target_se_intron_up.tsv",sep="\t",as.is=T,header=T,row.names=1);
clip_binding_data_intron_up[,"exonId"]<-str_sub(clip_binding_data_intron_up[,"exonId"],6,-1);

clip_names_intron_up<-paste0(colnames(clip_binding_data_intron_up)[6:ncol(clip_binding_data_intron_up)],"_intron_up");
colnames(clip_binding_data_intron_up)[6:ncol(clip_binding_data_intron_up)]<-clip_names_intron_up


clip_binding_data_intron_down<-read.table("../anno/gene_clip_target_se_intron_down.tsv",sep="\t",as.is=T,header=T,row.names=1);
clip_binding_data_intron_down[,"exonId"]<-str_sub(clip_binding_data_intron_down[,"exonId"],6,-1);

clip_names_intron_down<-paste0(colnames(clip_binding_data_intron_down)[6:ncol(clip_binding_data_intron_down)],"_intron_down");
colnames(clip_binding_data_intron_down)[6:ncol(clip_binding_data_intron_down)]<-clip_names_intron_down

clip_binding_data<-inner_join(clip_binding_data_exon,clip_binding_data_intron_down,by=c("gene_name"="gene_name","functions"="functions","exonId"="exonId"));
clip_binding_data<-inner_join(clip_binding_data,clip_binding_data_intron_up,by=c("gene_name"="gene_name","functions"="functions","exonId"="exonId"));

clip_names<-c(clip_names_intron_down,clip_names_intron_up,clip_names_exon);
#write.table(clip_binding_data[,c("gene_name","functions","exonId",clip_names) ],file="../anno/clip_names_binding_data.tsv",sep="\t",quote=F,row.names=F);


expression_psi_data<-read.table(file="../anno/all_data_expression_psi.tsv",as.is=T,header=T);
sample_names<-colnames(expression_psi_data);

expression_psi_data<-cbind(rownames(expression_psi_data),expression_psi_data);

colnames(expression_psi_data)[1]<-"name";
expression_psi_data[,"name"]<-as.character(expression_psi_data[,"name"]);

clip_binding_exp_psi_data<-inner_join(clip_binding_data,expression_psi_data,by=c("exonId"="name") );

#remove event which share same first and third exon, this will introduce bias.
#first_exon_start<-sapply(strsplit(clip_binding_exp_psi_data[,"exonId"],":"),"[",2 );
#first_exon_end<-sapply(strsplit(clip_binding_exp_psi_data[,"exonId"],":"),"[",3 );
#third_exon_start<-sapply(strsplit(clip_binding_exp_psi_data[,"exonId"],":"),"[",8 );
#third_exon_end<-sapply(strsplit(clip_binding_exp_psi_data[,"exonId"],":"),"[",9 );


non_duplicated_exon_ids<-!duplicated(cbind(first_exon_start,first_exon_end,third_exon_start,third_exon_end));

write.table(clip_binding_exp_psi_data,file="../anno/clip_binding_psi.tsv",sep="\t",row.names=F,col.names=T);

median_psi<-apply(clip_binding_exp_psi_data[,sample_names ],1,function(x){
  return(median(as.numeric(x),na.rm=T) );
});

var_psi<-apply(clip_binding_exp_psi_data[,sample_names ],1,function(x){
  return(var(as.numeric(x),na.rm=T) );
});

count_psi<-apply(clip_binding_exp_psi_data[,sample_names ],1,function(x){
  return(sum(!is.na(x) ) );
});

clip_binding_exp_psi_data<-cbind(clip_binding_exp_psi_data,median_psi,var_psi,count_psi);

clip_binding_median_psi_info<-clip_binding_exp_psi_data[,c("gene_name","exonId","median_psi","var_psi","count_psi",clip_names )];
write.table(clip_binding_median_psi_info,file="../anno/clip_binding_psi_data.tsv",sep="\t",quote=F,row.names=F);


clip_binding_median_psi<-clip_binding_exp_psi_data[,c(clip_names, "median_psi" )];

clip_binding_median_psi_lm<-lm(median_psi~.,data=clip_binding_median_psi);
clip_binding_median_psi_lm_sum<-summary(clip_binding_median_psi_lm);

clip_binding_median_psi_lm_pValue<--1*log10(clip_binding_median_psi_lm_sum$coefficients[,4]);
clip_binding_median_psi_lm_coe<-clip_binding_median_psi_lm_sum$coefficients[,1];


#clip_binding_median_psi_lm_pValue<-sort(clip_binding_median_psi_lm_pValue)[1:(length(clip_binding_median_psi_lm_pValue)-1)];

pdf("../anno/linear regression of RBP(exon_intron) and Event.pdf",width=8,height=5);
#barplot(clip_binding_median_psi_lm_pValue,las=2,main="Linear regression P-value of RBP and Event(-log10)");

#plot(clip_binding_median_psi_lm_coe,clip_binding_median_psi_lm_pValue);
#text(clip_binding_median_psi_lm_coe,clip_binding_median_psi_lm_pValue,labels=names(clip_binding_median_psi_lm_coe) );
p<-ggplot()+geom_point(aes(x=clip_binding_median_psi_lm_coe[-1],y=clip_binding_median_psi_lm_pValue[-1]),shape=1 ,size=0.01)+
geom_text(aes(x=clip_binding_median_psi_lm_coe[-1],y=clip_binding_median_psi_lm_pValue[-1],label=names(clip_binding_median_psi_lm_coe)[-1]),size=1 )+#scale_y_log10()
#scale_x_log10(breaks=c(0,1e-05,1e-04,1e-03,1e-02,1e+1,1))+scale_y_log10(),;
xlab("coefficients")+ylab("p-value(-log10)");

print(p);
dev.off();


clip_binding<-clip_binding_median_psi[,-ncol(clip_binding_median_psi)];
#clip_binding_cor_mat<-t(as.matrix(clip_binding));
clip_binding_cor_mat<-Matrix(as.matrix(clip_binding) );

clip_binding_cor<-corSparse(clip_binding_cor_mat);

colnames(clip_binding_cor)<-colnames(clip_binding)
rownames(clip_binding_cor)<-colnames(clip_binding)

pdf("../anno/RBP(exon_intron)_cor.pdf");
pheatmap(clip_binding_cor,cluster_rows=T,cluster_cols=T,fontsize_row=3,fontsize_col=3);
dev.off();



#clip_binding_cor_binary<-sapply(clip_binding_cor,function(x){if(x<0)return(-1);return(1); } );
#pdf("../anno/RBP(exon_intron)_cor_binary.pdf");
#pheatmap(clip_binding_cor_binary,cluster_rows=F,cluster_cols=F);
#dev.off();
