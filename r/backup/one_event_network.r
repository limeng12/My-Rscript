
clip_exon_names<-c(clip_names,
                   "chr11:8705553:8705628:+@chr11:8706233:8706439:+@chr11:8707225:8707395:+");

one_event_clip_data<-as.data.frame(t(all_data[clip_exon_names,]) ); 
one_event_clip_data_naRemove<-
  one_event_clip_data[!is.na(one_event_clip_data[,ncol(one_event_clip_data)]), ]; 

name<-colnames(one_event_clip_data_naRemove)[ncol(one_event_clip_data_naRemove)]; 
colnames(one_event_clip_data_naRemove)[ncol(one_event_clip_data_naRemove)]<-"event"; 


library(Rgraphviz); 
bn<-gs(one_event_clip_data_naRemove); 

library(deal); 
ne<-network(one_event_clip_data_naRemove); 
ne.prior<-jointprior(ne); 
rne<-learn(ne,one_event_clip_data_naRemove); 


pdf("../anno/bayesian_network.pdf"); 
plot(bn); 
graphviz.plot(bn); 
dev.off(); 
#
