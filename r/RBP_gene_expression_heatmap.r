setwd("E:\\limeng\\RBP_network\\r");

all_clip_event_data<-read.table(file="../anno/all_data_expression_psi.tsv",as.is=TRUE,header=TRUE);

all_data_label<-apply(all_clip_event_data,1,
                      function(x){if( sum(!is.na(x)) ==0 ) return(FALSE);return(TRUE);  } )

all_data<-all_clip_event_data[all_data_label,];
correMat<-matrix(ncol=48,nrow=(nrow(all_data)-48) );

for(i in 1:48){
  print(paste0("RNA binding factort: ", i) ); 
  
  correMat[,i]<-abs(cor(as.numeric(all_data[i,]),t(all_data[(49):nrow(all_data),]),
                    method="spearman",use="pairwise.complete.obs") ); 
  
}

rownames(correMat)<-rownames(all_data)[49:nrow(all_data)]; 
colnames(correMat)<-rownames(all_data)[1:48]; 

a<-apply(correMat,2,function(x){sum(x,na.rm=TRUE)})
labels<-names(sort(a) );

#correMat<-correMat[order(apply(correMat,1,sum) ), ]; 

#write.table(correMat,file="../anno/clip_event_correlation_matrix.tsv",sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE,); 

#gene_clip_target<-read.table("../anno/gene_clip_target_se_exon.tsv",header=TRUE,as.is=TRUE,sep="\t"); 
RBP_binding_psi_data<-read.table("../anno/clip_binding_medianPsi.tsv",header=TRUE,as.is=TRUE,sep="\t");
filter_sample_names<-readLines("../anno/filter_samples.tsv");

#RBP_binding_data<-RBP_binding_psi_data[,1:153];

######################################################filter intron and sample#############################################################################
fileter_psi_names<-colnames(RBP_binding_psi_data)[ grepl("EN",colnames(RBP_binding_psi_data) )];
#|intron
RBP_binding_data<-RBP_binding_psi_data[,setdiff(colnames(RBP_binding_psi_data),c("median_psi",filter_sample_names,fileter_psi_names) ) ];

RBP_binding_data_sum<-apply(RBP_binding_data[,2:ncol(RBP_binding_data)],1,function(x){sum(x!=0) });
gene_clip_target<-cbind(RBP_binding_data_sum,RBP_binding_data);

colnames(gene_clip_target)[1]<-"Sum";


#get the clip target and event name
library(stringr);
high_dense_clip_events<- gene_clip_target[gene_clip_target[,"Sum"]>70,"exonId"];
#high_dense_clip_events<-str_sub(high_dense_clip_events,6);

low_dense_clip_events<- gene_clip_target[gene_clip_target[,"Sum"]<5,"exonId"];
#low_dense_clip_events<-str_sub(low_dense_clip_events,6); 

correMat_high_clip<-correMat[intersect(high_dense_clip_events,rownames(correMat) ),]; 
correMat_low_clip<-correMat[intersect(low_dense_clip_events,rownames(correMat) ),]; 

#remove NAs
correMat_high_clip_label<-apply(correMat_high_clip,1,
                                function(x){if( sum(!is.na(x)) ==0 ) return(FALSE);return(TRUE);  } )
correMat_high_clip<-correMat_high_clip[correMat_high_clip_label,];

correMat_low_clip_label<-apply(correMat_low_clip,1,
                               function(x){if( sum(!is.na(x)) ==0 ) return(FALSE);return(TRUE);  } )
correMat_low_clip<-correMat_low_clip[correMat_low_clip_label,];

write.table(correMat_high_clip,file="../result/high_clip_corr_mat.tsv",sep="\t"); 
write.table(correMat_low_clip,file="../result/low_clip_corr_mat.tsv",sep="\t"); 


get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

#upper_tri_correMat_low_clip<-get_upper_tri(correMat_low_clip[1:nrow(correMat_high_clip),]);
#upper_tri_correMat_high_clip<-get_upper_tri(correMat_high_clip[,labels]);

library(reshape2)

correMat_low_clip_sample<-correMat_low_clip[sample(1:nrow(correMat_low_clip), nrow(correMat_high_clip) ), ];


source("multiplot.r")

library(gplots)
library(pheatmap); 
library(ggplot2)

melt_upper_tri_correMat_low_clip<-melt(correMat_low_clip);
melt_upper_tri_correMat_high_clip<-melt(correMat_high_clip);


#correMat_high_clip[1,1]<-1;correMat_high_clip[1,10]<-0;
#pheatmap(correMat_high_clip[,labels],cluster_rows=T,cluster_cols=F,show_rownames=F ,main="Events with RBP numbers > 40"); 
#pheatmap(correMat_low_clip[1:nrow(correMat_high_clip),labels],cluster_rows=T,cluster_cols=F ,show_rownames=F,main="Events with RBP numbers < 5" ); 

p1<-ggplot(data = melt_upper_tri_correMat_low_clip, aes(Var2, Var1, fill = value))+
 geom_tile()+
 scale_fill_gradient2(low = "white", high = "red", mid = "yellow", 
   midpoint = 0.5, limit = c(0,1), #space = "Lab", 
   name="Spearman\nCorrelation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 90, vjust = 1, 
    size = 6, hjust = 1),axis.text.y = element_blank())+ggtitle("Events' PSI correlate with RBPs' expression")+
 xlab("RBPs")+ylab("Events with number of RBP bindings < 5")+guides(fill=FALSE)

p2<-ggplot(data = melt_upper_tri_correMat_high_clip, aes(Var2, Var1, fill = value))+
 geom_tile()+
 scale_fill_gradient2(low = "white", high = "red", mid = "yellow", 
   midpoint = 0.5, limit = c(0,1),# space = "Lab", 
   name="Spearman\nCorrelation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 90, vjust = 1, 
    size = 6, hjust = 1),axis.text.y = element_blank() )+ggtitle("")+
 xlab("RBPs")+ylab("Events with number of RBP bindings > 40") + theme()#+guide_legend(title="Correlation\ncoefficient")#+guides(fill=FALSE);


#multiplot(p1,p2,cols=2,width=c(0.78,1));

jpeg("../result/Figure 2__RBP_event_correlation_matrix1.jpeg",width=1200,height=1200);
pheatmap(correMat_high_clip,cluster_rows=TRUE,cluster_cols=FALSE,show_rownames=FALSE,fontsize_col=3 ,main="Events with RBP numbers > 30");
dev.off();


jpeg("../result/Figure 2__RBP_event_correlation_matrix2.jpeg",width=1200,height=1200);
pheatmap(correMat_low_clip_sample,cluster_rows=TRUE,cluster_cols=FALSE ,show_rownames=FALSE,fontsize_col=3,main="Events with RBP numbers < 5" );
dev.off(); 



