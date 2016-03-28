setwd("/Users/mengli/Documents/projects/RBP_network/r");
library(pheatmap)
library(ggplot2)
library(dplyr)
pdf("../result/qc.pdf",width=10,height=6)

#postscript("../result/qc.ps")
RBP_binding_psi_data<-read.table("../anno/clip_binding_medianPsi.tsv",header=TRUE,as.is=TRUE,sep="\t");
#filter_sample_names<-readLines("../anno/filter_samples.tsv");
#RBP_binding_data<-RBP_binding_psi_data[,1:153];
#fileter_psi_names<-colnames(RBP_binding_psi_data)[ grepl("EN",colnames(RBP_binding_psi_data) )];
#RBP_binding_data<-RBP_binding_psi_data[,setdiff(colnames(RBP_binding_psi_data),c("start.y","end.y","start","end",filter_sample_names,fileter_psi_names) ) ];
RBP_binding_data<-RBP_binding_psi_data[,-2];

RBP_binding_sites_sum<-apply(RBP_binding_data[,1:ncol(RBP_binding_data)],1,sum);
names(RBP_binding_sites_sum)<-colnames(RBP_binding_data)[1:ncol(RBP_binding_data)];


#hist(RBP_binding_sites_sum,main="Number of RBPs binding on each event",breaks=50);
RBP_binding_sites_sum_plot<-ggplot() + geom_histogram(aes(x=RBP_binding_sites_sum ),bins=50 )+
scale_x_log10()+
ggtitle("Number of RBPs binding sites in each event") ;
print(RBP_binding_sites_sum_plot);
library(pheatmap)
########
setwd("/Users/mengli/Documents/projects/RBP_network/r");
library(pheatmap)

gene_expression_data<-read.table("../anno/gene_expression_table.tsv",header=TRUE,as.is=TRUE);
miso_psi_data<-read.table("../anno/miso_event.tsv",header=TRUE,as.is=TRUE);
filter_sample_names<-readLines("../anno/filter_samples.tsv");

overlap_sample_names<-intersect(colnames(miso_psi_data),colnames(gene_expression_data) );
overlap_sample_names<-setdiff(overlap_sample_names,filter_sample_names);

#gene_expression_names_correct<-c(overlap_sample_names);
#miso_psi_names_correct<-c(overlap_sample_names);


gene_expression_data<-gene_expression_data[!duplicated(gene_expression_data[,"X.gene"]),];
rownames(gene_expression_data)<-gene_expression_data[,"X.gene"];
rownames(miso_psi_data)<-miso_psi_data[,"event_name"];

gene_expression_data<-gene_expression_data[,overlap_sample_names];
miso_psi_data<-miso_psi_data[,overlap_sample_names];

#miso_psi_data<-miso_psi_data[,-ncol(miso_psi_data)];
miso_psi_data_melt_sampleAvg<-apply(miso_psi_data,2,function(x){median(x,na.rm=TRUE) } );
miso_psi_data_sampleSum<-apply(miso_psi_data,2,function(x){ sum(!is.na(x) ) } );
names(miso_psi_data_sampleSum)<-colnames(miso_psi_data)

miso_psi_data_sampleSum<-sort(miso_psi_data_sampleSum)
barplot(miso_psi_data_sampleSum,las=2,main="Each's sample's Events number",cex.names=0.5);


library(reshape2)
miso_psi_data_melt<-melt(miso_psi_data);

miso_psi_data_melt_noNa<-miso_psi_data_melt[!is.na(miso_psi_data_melt[,"value"]),];

miso_psi_data_melt_noNa[,"variable"]<-factor(miso_psi_data_melt_noNa[,"variable"],
                                             levels=colnames(miso_psi_data)[order(miso_psi_data_melt_sampleAvg) ] )

miso_psi_sampleSum_data<-as.data.frame(table(miso_psi_data_melt_noNa[,"variable"]) );
############################################################too many missing values, make heatmap useless##############################################
#miso_psi_data<-as.matrix(miso_psi_data[sample(nrow(miso_psi_data),10000),]);
#var_1<-(apply(miso_psi_data,1,function(x){return(var(x,na.rm=TRUE)) } )<0.1 )
#var_2<-(apply(miso_psi_data,2,function(x){return(var(x,na.rm=TRUE)) } )<0.1 )
#miso_psi_data<-miso_psi_data[var_1,];
#var_2[length(var_2)]<-TRUE;
#miso_psi_data<-miso_psi_data[,var_2];
#miso_psi_data_naLabel_1<-apply(miso_psi_data,1,function(x){return(sum(!is.na(x))>3) } );
#miso_psi_data<-miso_psi_data[miso_psi_data_naLabel_1,];
#miso_psi_data_naLabel_2<-apply(miso_psi_data,2,function(x){return(sum(!is.na(x))>3) } );
#miso_psi_data<-miso_psi_data[,miso_psi_data_naLabel_2];
########################################################################################################################################################
miso_psi_data_melt_noNa_join<-inner_join(miso_psi_data_melt_noNa,miso_psi_sampleSum_data,by=c("variable"="Var1") );

miso_pis_box_plot<-ggplot(miso_psi_data_melt_noNa_join,aes(x=variable, y=value  ))+geom_boxplot(aes(fill=Freq)  ) +theme_classic()+
scale_fill_gradient(trans = "log",low="white",high="red")+
ylab("Miso Events' PSI")+
theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 6, hjust = 1),axis.text.y = element_blank() )
print(miso_pis_box_plot)

#gene_expression_data<-gene_expression_data[,-1];
gene_expression_data_sampleAvg<-apply(gene_expression_data,2,function(x){median(log10(x),na.rm=TRUE) } );
gene_expression_data_sampleSum<-apply(gene_expression_data,2,function(x){ sum(!is.na(x) ) } );
names(gene_expression_data_sampleSum)<-colnames(gene_expression_data);

gene_expression_data_sampleSum<-sort(gene_expression_data_sampleSum);
barplot(gene_expression_data_sampleSum,las=2,main="Each's sample's Gene expression number",cex.names=0.5);

gene_expression_data_melt<-melt(gene_expression_data);
gene_expression_data_melt[,"variable"]<-factor(gene_expression_data_melt[,"variable"],
                                               levels=colnames(gene_expression_data)[order(gene_expression_data_sampleAvg) ] );
gene_expression_data_melt_noNa<-gene_expression_data_melt[!is.na(gene_expression_data_melt[,"value"]),];

gene_exp_box_plot<-ggplot(gene_expression_data_melt_noNa)+
  geom_boxplot(aes(x=variable, y=value  ) ) +
  theme_classic()+
  ylab("Gene expression' FPKM")+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 6, hjust = 1),
        axis.text.y = element_text( vjust = 1, size = 6, hjust = 1) )+
  scale_y_log10();

print(gene_exp_box_plot);


ggene_expression_data_sample<-gene_expression_data[sample(nrow(gene_expression_data),5000 ) ,-1];

ggene_expression_data_sample<-log(ggene_expression_data_sample+0.001)


pheatmap(ggene_expression_data_sample,
  cluster_rows=TRUE,cluster_cols=TRUE,show_rownames=FALSE,show_colnames=TRUE,
  fontsize_col=5 ,main="Gene expression of 1000 genes heatmap (log)");

dev.off();

if(FALSE){
  jpeg("../result/gene_expression_box_plot.jpeg",width=1000,height=800);
  print(gene_exp_box_plot);
  dev.off();
  
}

