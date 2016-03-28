setwd("E:\\limeng\\RBP_network\\r");

clip_binding_psi_data<-read.table(file="../anno/clip_binding_psi_data.tsv",sep="\t",header=TRUE,as.is=TRUE);


clip_binding_data<-clip_binding_psi_data[,6:ncol(clip_binding_psi_data)];

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


colnames(RBP_bind_region_specific)<-c("RBP bind intron upstream","RBP bind intron downstream","RBP bind exon")
rownames(RBP_bind_region_specific)<-RBP_names

library(reshape2)
RBP_bind_region_specific_melt<-melt(as.matrix(RBP_bind_region_specific,na.rm=T) );

RBP_bind_region_specific_melt[,"value"]<-log10(RBP_bind_region_specific_melt[,"value"])

p1<-ggplot(data = RBP_bind_region_specific_melt, aes(Var2, Var1, fill = value) )+
 geom_tile()+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint =  summary(RBP_bind_region_specific_melt[,"value"])[2],
   limit = c(min(RBP_bind_region_specific_melt[,"value"]),max(RBP_bind_region_specific_melt[,"value"])), #space = "Lab", 
   name="Number of \nbinding sites\n(log10)") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 0, vjust = 1, 
    size = 10, hjust = 0.5),axis.text.y = element_text(angle = 0, vjust = 0.6, 
    size = 10, hjust = 0))+
 xlab("")+ylab("RBPs")+ggtitle("RRBP binding region specific")#+guides(fill=TRUE)

print(p2)


