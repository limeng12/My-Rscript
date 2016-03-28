setwd("E:\\limeng\\RBP_network\\r");

clip_binding_psi_data<-read.table(file="../anno/clip_binding_psi_var_median.tsv",sep="\t",header=TRUE,as.is=TRUE);

#filter
clip_binding_psi_data_filter<-subset(clip_binding_psi_data, (median_psi>0.9|median_psi<0.1)&var_psi<0.1&count_psi>70 );

is_constitute<-sapply(clip_binding_psi_data_filter[,"median_psi"],function(x){if(x>0.9) return(TRUE);return(FALSE) } )

clip_binding_psi_data_filter<-cbind(is_constitute,clip_binding_psi_data_filter);

clip_binding_psi_data_filter_class<-clip_binding_psi_data_filter[,-c(2:6)];

cor_of_RBP_psi<-cor(clip_binding_psi_data_filter_class,is_constitute,method="spearman");

cor_of_RBP_names<-rownames(cor_of_RBP_psi);
intron_up_names<-cor_of_RBP_names[ grepl("up",cor_of_RBP_names)];
intron_down_names<-cor_of_RBP_names[ grepl("down",cor_of_RBP_names)];
exon_names<-cor_of_RBP_names[ grepl("exon",cor_of_RBP_names)];

RBP_binding_pos_trend<-data.frame(bind_intron_upStream=cor_of_RBP_psi[intron_up_names,],
bind_intron_downStream=cor_of_RBP_psi[intron_down_names,],
bind_exon=cor_of_RBP_psi[exon_names,]
);


RBP_names<-sapply(strsplit(exon_names,"_"),"[",1);
rownames(RBP_binding_pos_trend)<-RBP_names

RBP_binding_pos_trend_mat<-as.matrix(RBP_binding_pos_trend);

colnames(RBP_binding_pos_trend_mat)<-c("intron upstream","intron downstream","exon")

#RBP_binding_pos_trend_mat<-apply(RBP_binding_pos_trend_mat,c(1,2),function(x){if(x>0)return(1);return(0) })

write.table(RBP_binding_pos_trend_mat,file="../result/RBPs' Binding region specifity.tsv",sep="\t",quote=FALSE)

library(gplots)
library(pheatmap);
library(grid)     ## Need to attach (and not just load) grid package
library(pheatmap)
## Edit body of pheatmap:::draw_colnames, customizing it to your liking
draw_colnames_45 <- function (coln, ...) {
    m = length(coln)
    x = (1:m)/m - 1/2/m
    grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
        hjust = 1, rot = 0, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
ns=asNamespace("pheatmap"))

#pheatmap(d)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
RBP_binding_pos_trend_mat_melt<-melt(RBP_binding_pos_trend_mat,na.rm=T);

#pheatmap(RBP_binding_pos_trend_mat,cluster_rows=F,cluster_cols=F,fontsize_row=5,
# fontsize_col=7,
# color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10),
 #breaks = c(-0.4,-0.3,-0.2,-0.1, 0 , 0.02, 0.05, 0.07, 0.1)
# breaks = c(-0.15, -0.105, -0.075, -0.035, 0 , 0.1, 0.2, 0.3, 0.43),
# main="correlation between number of RBPs and PSI"
# );
p1<-ggplot(data = RBP_binding_pos_trend_mat_melt, aes(Var2, Var1, fill = value) )+
 geom_tile()+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(min(RBP_binding_pos_trend_mat_melt[,"value"]),max(RBP_binding_pos_trend_mat_melt[,"value"])), #space = "Lab", 
   name="Spearman\nCorrelation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 0, vjust = 1, 
    size = 6, hjust = 0.5),axis.text.y = element_text(angle = 0, vjust = 0.6, 
    size = 10, hjust = 0))+
 xlab("")+ylab("RBPs")+ggtitle("(A) RBP-region combinations correlate with PSI")#+guides(fill=TRUE)


setwd("E:\\limeng\\RBP_network\\r");

clip_binding_psi_data<-read.table(file="../anno/clip_binding_psi_var_median.tsv",sep="\t",header=TRUE,as.is=TRUE);


clip_binding_data<-clip_binding_psi_data[,-c(1:5)];

clip_binding_data_sum<-apply(clip_binding_data,2,sum);

cor_of_RBP_names<-names(clip_binding_data_sum);
intron_up_names<-cor_of_RBP_names[ grepl("up",cor_of_RBP_names)];
intron_down_names<-cor_of_RBP_names[ grepl("down",cor_of_RBP_names)];
exon_names<-cor_of_RBP_names[ grepl("exon",cor_of_RBP_names)];


RBP_names<-sapply(strsplit(exon_names,"_"),"[",1);

RBP_bind_region_specific<-data.frame(bind_intron_upStream=clip_binding_data_sum[intron_up_names],
bind_intron_downStream=clip_binding_data_sum[intron_down_names],
bind_exon=clip_binding_data_sum[exon_names]
);

colnames(RBP_bind_region_specific)<-c("intron upstream","intron downstream","exon")
rownames(RBP_bind_region_specific)<-RBP_names

library(reshape2)
RBP_bind_region_specific_melt<-melt(as.matrix(RBP_bind_region_specific,na.rm=T) );

RBP_bind_region_specific_melt[,"value"]<-log10(RBP_bind_region_specific_melt[,"value"])

p2<-ggplot(data = RBP_bind_region_specific_melt, aes(Var2, Var1, fill = value) )+
 geom_tile()+
 scale_fill_gradient2(low = "white", high = "red", mid = "yellow", midpoint =  summary(RBP_bind_region_specific_melt[,"value"])[2],
   limit = c(min(RBP_bind_region_specific_melt[,"value"]),max(RBP_bind_region_specific_melt[,"value"])), #space = "Lab", 
   name="Number of \nbinding sites\n(log10)") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 0, vjust = 1, 
    size = 6, hjust = 0.5),axis.text.y = element_text(angle = 0, vjust = 0.6, 
    size = 10, hjust = 0))+
 xlab("")+ylab("RBPs")+ggtitle("(B) RBP binding region specific")#+guides(fill=TRUE)


source("multiplot.r");
pdf("../result/Figure 5__RBP_binding_pos_trend.pdf",width=10,height=7);
multiplot(p1,p2,cols=2);

dev.off();
