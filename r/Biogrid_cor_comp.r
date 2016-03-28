library(readr)
setwd("E:\\limeng\\RBP_network\\r")

biogrid_data<-read_tsv("../anno/BIOGRID-ALL-3.4.134.tab.txt",skip=35,col_names=TRUE);
clip_name_map<-read.table("../anno/clipNameTargetMap.txt",sep="\t",header=FALSE,as.is=TRUE);
biogrid_data_data<-as.data.frame(biogrid_data,stringsAsFactors=FALSE);

rbp_names<-unique(clip_name_map[,1])

gene_interact<-biogrid_data_data[,c("OFFICIAL_SYMBOL_A","OFFICIAL_SYMBOL_B")];

rbp_interact<-matrix(ncol=2,nrow=0);

for(i in 1:nrow(gene_interact) ){
  gene_symbol1<-gene_interact[i,1];
  gene_symbol2<-gene_interact[i,2];
  if(gene_symbol1==gene_symbol2)
    next;
  
  if(is.element(gene_symbol1,rbp_names)&&is.element(gene_symbol2,rbp_names) ){
   print( paste0(gene_symbol1,"\t",gene_symbol2) );
   rbp_interact<-rbind(rbp_interact,c(gene_symbol1,gene_symbol2) );
  }
}

rbp_interact_u<-unique(rbp_interact);
rbp_interact2<-rbp_interact_u;


rbp_interact2[,1]<-rbp_interact_u[,2]
rbp_interact2[,2]<-rbp_interact_u[,1]

rbp_interact<-rbind(rbp_interact_u,rbp_interact2);
rbp_interact_unique<-unique(rbp_interact);

rbp_interact_unique<-cbind(rbp_interact_unique,rep(TRUE,nrow(rbp_interact_unique)) );
colnames(rbp_interact_unique)<-c("Var1","Var2","T");

rbp_cor<-read.table(file="../result/RBP_exon_binding_cor_spearman.tsv",sep="\t",header=TRUE,as.is=TRUE );
rbp_cor_data<-rbp_cor[!(rbp_cor[,"Var1"]==rbp_cor[,"Var2"]),];

rbp_cor_data_sort<-rbp_cor_data[order(rbp_cor_data[,"value"],decreasing=TRUE),]
rbp_cor_data_sort<-unique(rbp_cor_data_sort);

name1<-sapply(strsplit(rbp_cor_data_sort[,"Var1"],"_"),"[",1)
name2<-sapply(strsplit(rbp_cor_data_sort[,"Var2"],"_"),"[",1)

rbp_cor_data_sort<-cbind(rbp_cor_data_sort,name1,name2);
colnames(rbp_cor_data_sort)<-c("RBP_region1","RBP_region2","Cor_effiients","Var1","Var2");

#rbp_cor_data_sort[1,];rbp_cor_data_sort[,2];

rbp_cor_data_biogrid<-left_join(rbp_cor_data_sort,rbp_interact_unique,by=c("Var1"="Var1","Var2"="Var2"),copy=TRUE );
rbp_cor_data_biogrid_u<-rbp_cor_data_biogrid[seq(1,nrow(rbp_cor_data_biogrid),by=2),];
write.table(rbp_cor_data_biogrid_u,file="../result/rbp_cor_biogrid.tsv",sep="\t",row.names=FALSE);


############################################################################################################################
library(igraph)

rbp_interact_unique_names<-apply(rbp_interact_unique,1,function(x){paste0(sort(x),collapse="-")})

rbp_interact_unique_one<-rbp_interact_unique[!duplicated(rbp_interact_unique_names),];
ne<-graph_from_edgelist(rbp_interact_unique_one[,1:2],directed=FALSE);

pdf("../anno/bio_grid_RBP_net.pdf")
plot(ne,layout=layout_as_star,vertex.size=3,vertex.label.cex=0.5);
dev.off();

ne_degree<-degree(ne);ne_degree_sort<-sort(ne_degree);
write.table(ne_degree_sort,file="../result/degree of freedom.tsv",sep="\t",row.names=TRUE )

