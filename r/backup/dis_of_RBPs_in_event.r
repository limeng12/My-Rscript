setwd("/home/limeng/RBP_network/r")
RBP_binding_psi_data<-read.table("../anno/clip_binding_psi.tsv",header=T,as.is=T);

RBP_binding_data<-RBP_binding_psi_data[,1:153];
RBP_binding_data<-RBP_binding_data[,setdiff(colnames(RBP_binding_data),c("start.y","end.y","start","end") ) ];

RBP_binding_sites_sum<-apply(RBP_binding_data[,6:ncol(RBP_binding_data)],1,sum);

