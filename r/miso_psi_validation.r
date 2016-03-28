setwd("E:\\limeng\\RBP_network\\r")
target_events<-read.table("../anno/RBP_binding_psi_for_logistic_regression.tsv",header=TRUE,sep="\t",as.is=TRUE);
filter_sample_names<-readLines("../anno/filter_samples.tsv");

target_events_exonIds<-target_events[,"exonId"];

miso_event_0count_data<-read.table("../anno/miso_event_0count.tsv",header=TRUE,sep="\t",as.is=TRUE);
miso_event_1count_data<-read.table("../anno/miso_event_1count.tsv",header=TRUE,sep="\t",as.is=TRUE)
miso_event_ci<-read.table("../anno/miso_event_ci.tsv",header=TRUE,sep="\t",as.is=TRUE)

miso_event_0count_data<-miso_event_0count_data[,setdiff(colnames(miso_event_0count_data),filter_sample_names)];
miso_event_1count_data<-miso_event_1count_data[,setdiff(colnames(miso_event_1count_data),filter_sample_names)];
miso_event_ci<-miso_event_ci[,setdiff(colnames(miso_event_ci),filter_sample_names)];


miso_event_0count_data_targetEvents<-miso_event_0count_data[is.element(miso_event_0count_data[,"event_name"],target_events_exonIds),];
miso_event_1count_data_targetEvents<-miso_event_1count_data[is.element(miso_event_1count_data[,"event_name"],target_events_exonIds),]
miso_event_ci_targetEvents<-miso_event_ci[is.element(miso_event_ci[,"event_name"],target_events_exonIds),]

miso_event_0count_data_targetEvents_avg<-apply(miso_event_0count_data_targetEvents[,-ncol(miso_event_0count_data_targetEvents)],1,function(x){median(x,na.rm=TRUE)} );
miso_event_1count_data_targetEvents_avg<-apply(miso_event_1count_data_targetEvents[,-ncol(miso_event_1count_data_targetEvents)],1,function(x){median(x,na.rm=TRUE)} );
miso_event_ci_targetEvents_avg<-apply(miso_event_ci_targetEvents[,-ncol(miso_event_ci_targetEvents)],1,function(x){median(x,na.rm=TRUE)} );

miso_event_counts<-data.frame(isoform0=miso_event_0count_data_targetEvents_avg,isoform1=miso_event_1count_data_targetEvents_avg);


source("multiplot.r")

pdf("../result/Supplement Figure 1__Miso_event_quality_control.pdf",width=10,height=6);
library(reshape2)
miso_event_counts_melt<-melt(miso_event_counts);
colnames(miso_event_counts_melt)<-c("Isoforms","reads counts")
#geom_boxplot
miso_event_counts_box_plot<-ggplot(miso_event_counts_melt)+geom_boxplot(aes(x=Isoforms, y=miso_event_counts_melt[,"reads counts"]  )  ) +
scale_y_log10()+theme_classic()+
ylab("Events' median reads count (log10)");

#print(miso_event_counts_box_plot);
#miso_event_0count_data_plot<-ggplot() +
#  geom_histogram(aes(x=miso_event_0count_data_targetEvents_avg) )+scale_x_log10();
  
#miso_event_1count_data_plot<-ggplot() +
#  geom_histogram(aes(x=miso_event_1count_data_targetEvents_avg) )+scale_x_log10();

median_ci<-median(miso_event_ci_targetEvents_avg);

miso_event_ci_data_plot<-ggplot() +
  geom_histogram(aes(x=miso_event_ci_targetEvents_avg,y=..density..),colour="black", fill="white" )+geom_density(color="blue")+
  ylab("frequency")+xlab("Events' median cofidence interval")+ggtitle(paste0("Median: ",median_ci ))+theme_classic();
  
multiplot(miso_event_counts_box_plot,miso_event_ci_data_plot,cols=2);

dev.off();

