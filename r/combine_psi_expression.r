setwd("/Users/mengli/Documents/projects/RBP_network/r");
library(pheatmap)
########
gene_expression_data<-read.table("../anno/gene_expression_table.tsv",header=TRUE,as.is=TRUE);
miso_psi_data<-read.table("../anno/miso_event.tsv",header=TRUE,as.is=TRUE);

filter_sample_names<-readLines("../anno/filter_samples.tsv");

overlap_sample_names<-intersect(colnames(miso_psi_data),colnames(gene_expression_data) );
overlap_sample_names<-setdiff(overlap_sample_names,filter_samples)

gene_expression_names_correct<-c(overlap_sample_names);
miso_psi_names_correct<-c(overlap_sample_names);

rownames(gene_expression_data)<-gene_expression_data[,"X.gene"];
rownames(miso_psi_data)<-miso_psi_data[,"event_name"];

#pdf("../result/heatmap of PSI and Expression")
#pheatmap(log10(gene_expression_data[sample(nrow(gene_expression_data),1000 ) ,-1]+0.1 ) ,
#cluster_rows=FALSE,cluster_cols=TRUE,show_rownames=FALSE,show_colnames=FALSE );

#pheatmap(miso_psi_data[sample(nrow(miso_psi_data),1000),-115],
#cluster_rows=FALSE,cluster_cols=TRUE,show_rownames=FALSE,show_colnames=FALSE);
#dev.off();

clip_names<-read.table("../anno/clipNameTargetMap.txt",header=FALSE,sep="\t",as.is=TRUE)[,1];
clip_names<-unique(clip_names);
#clip_exon_names<-c(clip_names,"chr3:9854932:9855029:+@chr3:9859329:9862425:+@chr3:9867484:9867632:+");


all_clip_event_data<-rbind(gene_expression_data[clip_names,gene_expression_names_correct],
                miso_psi_data[,miso_psi_names_correct] );

all_data<-all_clip_event_data;

#write.table(all_data,file="../anno/all_data_expression_psi.tsv",
#            sep="\t",row.names = TRUE,col.names = TRUE,quote=FALSE); 

write.table(all_data,file="../result/all_data_expression_psi.tsv",
            sep="\t",row.names = TRUE,col.names = TRUE,quote=FALSE); 
            
#one_event<-all_data[clip_exon_names,];
#write.table(one_event,file="../anno/one_event_clip.tsv",
#sep="\t",row.names = TRUE, col.names = TRUE,quote = FALSE); 

