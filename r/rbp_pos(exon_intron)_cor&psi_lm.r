library(dplyr)
library(stringr)
library(Matrix)
#library(qlcMatrix)
library(pheatmap)
library(ppcor)
library(reshape2)
library(grid)     ## Need to attach (and not just load) grid package
library(pheatmap)
library(ggplot2)

setwd("/Users/mengli/Documents/projects/RBP_network/r");

clip_binding_data_exon<-read.table("../anno/gene_clip_target_se_exon.tsv",sep="\t",as.is=TRUE,header=TRUE,row.names=1);
clip_binding_data_exon[,"exonId"]<-str_sub(clip_binding_data_exon[,"exonId"],6,-1);

clip_names_exon<-paste0(colnames(clip_binding_data_exon)[6:ncol(clip_binding_data_exon)],"_exon");
colnames(clip_binding_data_exon)[6:ncol(clip_binding_data_exon)]<-clip_names_exon


clip_binding_data_intron_up<-read.table("../anno/gene_clip_target_se_intron_up.tsv",sep="\t",as.is=TRUE,header=TRUE,row.names=1);
clip_binding_data_intron_up[,"exonId"]<-str_sub(clip_binding_data_intron_up[,"exonId"],6,-1);

clip_names_intron_up<-paste0(colnames(clip_binding_data_intron_up)[6:ncol(clip_binding_data_intron_up)],"_intron_up");
colnames(clip_binding_data_intron_up)[6:ncol(clip_binding_data_intron_up)]<-clip_names_intron_up


clip_binding_data_intron_down<-read.table("../anno/gene_clip_target_se_intron_down.tsv",sep="\t",as.is=TRUE,header=TRUE,row.names=1);
clip_binding_data_intron_down[,"exonId"]<-str_sub(clip_binding_data_intron_down[,"exonId"],6,-1);

clip_names_intron_down<-paste0(colnames(clip_binding_data_intron_down)[6:ncol(clip_binding_data_intron_down)],"_intron_down");
colnames(clip_binding_data_intron_down)[6:ncol(clip_binding_data_intron_down)]<-clip_names_intron_down

clip_binding_data<-inner_join(clip_binding_data_exon,clip_binding_data_intron_down,by=c("gene_name"="gene_name","functions"="functions","exonId"="exonId"));
clip_binding_data<-inner_join(clip_binding_data,clip_binding_data_intron_up,by=c("gene_name"="gene_name","functions"="functions","exonId"="exonId"));

clip_names<-c(clip_names_intron_down,clip_names_intron_up,clip_names_exon);
#write.table(clip_binding_data[,c("gene_name","functions","exonId",clip_names) ],file="../anno/clip_names_binding_data.tsv",sep="\t",quote=FALSE,row.names=FALSE);
cat(clip_names,file="../anno/clip_names.txt",sep="\n")

expression_psi_data<-read.table(file="../anno/all_data_expression_psi.tsv",as.is=TRUE,header=TRUE);
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


#non_duplicated_exon_ids<-!duplicated(cbind(first_exon_start,first_exon_end,third_exon_start,third_exon_end));

#write.table(clip_binding_exp_psi_data,file="../anno/clip_binding_psi.tsv",sep="\t",row.names=FALSE,col.names=TRUE);

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
write.table(clip_binding_median_psi_info,file="../anno/clip_binding_psi_var_median.tsv",sep="\t",quote=FALSE,row.names=FALSE);

clip_binding_median_psi<-clip_binding_exp_psi_data[,c("median_psi",clip_names )];

clip_binding_median_psi_exon<-clip_binding_exp_psi_data[,c("median_psi","exonId",clip_names )];
write.table(clip_binding_median_psi_exon,file="../anno/clip_binding_medianPsi.tsv",sep="\t",row.names=FALSE,col.names=TRUE);

clip_binding_median_psi_event<-clip_binding_exp_psi_data[,c("exonId","median_psi",clip_names )];


event_split<-strsplit(clip_binding_median_psi_event[,"exonId"],":");
event_len<-c();
for(i in 1:length(event_split) ){
	if(event_split[[i]][10]=="+"){
		event_len[i]<-as.numeric(event_split[[i]][9])-as.numeric(event_split[[i]][2] );

	}else{
		event_len[i]<-as.numeric(event_split[[i]][3])-as.numeric(event_split[[i]][8] );

	}
}


clip_binding_median_psi_event<-cbind(event_len,clip_binding_median_psi_event);

##############################################################################Linear regression#############################################
clip_binding_median_psi<-clip_binding_exp_psi_data[,c("median_psi",clip_names )];
clip_binding_median_psi_lm<-lm(median_psi~.,data=clip_binding_median_psi);
clip_binding_median_psi_lm_sum<-summary(clip_binding_median_psi_lm);

pdf("../result/Figure 3__Histogram of residues.pdf",width=8,height=6)
hist(clip_binding_median_psi_lm$residuals,breaks=30, main="histogram of residues", xlab="residues");
dev.off();


#clip_binding_median_psi_lm_pValue<--1*log10(clip_binding_median_psi_lm_sum$coefficients[,4]);
#clip_binding_median_psi_lm_coe<-clip_binding_median_psi_lm_sum$coefficients[,1];

coe_logit<-clip_binding_median_psi_lm_sum$coefficients[order(clip_binding_median_psi_lm_sum$coefficients[,4]),]

logit_coeffcients<-coe_logit[2:21,1]
logit_coeffcients<-signif(logit_coeffcients,digits=4)
logit_p_value<-coe_logit[2:21,4]
logit_p_value<--log10(logit_p_value);
#clip_binding_median_psi_lm_pValue<-sort(clip_binding_median_psi_lm_pValue)[1:(length(clip_binding_median_psi_lm_pValue)-1)];

#pdf("../anno/linear regression of RBP(exon_intron) and Event.pdf",width=8,height=5);
#barplot(clip_binding_median_psi_lm_pValue,las=2,main="Linear regression P-value of RBP and Event(-log10)");
#plot(clip_binding_median_psi_lm_coe,clip_binding_median_psi_lm_pValue);
#text(clip_binding_median_psi_lm_coe,clip_binding_median_psi_lm_pValue,labels=names(clip_binding_median_psi_lm_coe) );
#p<-ggplot()+geom_point(aes(x=clip_binding_median_psi_lm_coe[-1],y=clip_binding_median_psi_lm_pValue[-1]),shape=1 ,size=0.01)+
#geom_text(aes(x=clip_binding_median_psi_lm_coe[-1],y=clip_binding_median_psi_lm_pValue[-1],label=names(clip_binding_median_psi_lm_coe)[-1]),size=1 )+#scale_y_log10()
#scale_x_log10(breaks=c(0,1e-05,1e-04,1e-03,1e-02,1e+1,1))+scale_y_log10(),;
#xlab("coefficients")+ylab("p-value(-log10)");
#print(p);

#logit_p_value_format_names<-names(logit_p_value_format)
#logit_p_value_format_names<-factor(logit_p_value_format_names,levels=logit_p_value_format_names );

#logit_p_value_format<-signif(logit_p_value,digits=3)
#p<-ggplot()+geom_bar(aes(x=logit_p_value_format_names,y=logit_p_value_format),stat="identity")+
# geom_text(aes(x=logit_p_value_format_names,y=logit_p_value_format,label=logit_coeffcients),vjust=-0.2,size=2)+theme_classic()+
# theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+ylab("p-value(-log10)")+ggtitle("top 20 P-values and coefficients for Linear regression")+
# xlab("RBP")
 #ylim(0,150)
#print(p);

#dev.off();


clip_binding<-clip_binding_median_psi[,clip_names];
#clip_binding_cor_mat<-t(as.matrix(clip_binding));
#clip_binding_cor_mat<-Matrix(as.matrix(clip_binding) );

#clip_binding_cor<-corSparse(clip_binding_cor_mat);
clip_binding_cor_spearman<-cor( sample_n(clip_binding,5000),method="spearman")

colnames(clip_binding_cor_spearman)<-colnames(clip_binding)
rownames(clip_binding_cor_spearman)<-colnames(clip_binding)

clip_binding_cor_pearson<-cor( sample_n(clip_binding,5000),method="pearson")

colnames(clip_binding_cor_pearson)<-colnames(clip_binding)
rownames(clip_binding_cor_pearson)<-colnames(clip_binding)

###############################partial correlation######################################################################
#partial_cor<-data.frame(Var1=c(),Var2=c(),value=c() );
partial_cor<-matrix(ncol=3,nrow=(ncol(clip_binding_median_psi_event)-3)*(ncol(clip_binding_median_psi_event)-3) )
#clip_binding_median_psi_event[1,(ncol(clip_binding_median_psi_event)-2)*(ncol(clip_binding_median_psi_event)-2)  ]
n<-1
for(i in 4:(ncol(clip_binding_median_psi_event)) ){
	print(i)
	for(j in 4:(ncol(clip_binding_median_psi_event)) ){
		var1<-colnames(clip_binding_median_psi_event)[i];
		var2<-colnames(clip_binding_median_psi_event)[j];
		
		tryCatch({
			lineIds<-sample(1:nrow(clip_binding_median_psi_event),5000);
			if(i!=j){
				est<-pcor.test(clip_binding_median_psi_event[1:lineIds,i],clip_binding_median_psi_event[1:lineIds,j],clip_binding_median_psi_event[1:lineIds,"event_len"],method="spearman")$estimate
			}else{
				est<--1;
			}
		}, error = function(e) {
			print(paste0(e," error here: ",var1," ",var2) );
		})
		partial_cor[n,]<-c(var1,var2,est);
		n<-n+1		
	}
	#if(i>2)
	#break;
	
}
colnames(partial_cor)<-c("Var1","Var2","value")
partial_cor<-as.data.frame(partial_cor)

partial_cor[,"Var1"]<-factor(partial_cor[,"Var1"],levels=colnames(clip_binding_median_psi_event)[4:ncol(clip_binding_median_psi_event) ] );
partial_cor[,"Var2"]<-factor(partial_cor[,"Var2"],levels=colnames(clip_binding_median_psi_event)[4:ncol(clip_binding_median_psi_event) ] );

partial_cor[,"value"]<-as.numeric(partial_cor[,"value"]);

max_cor<-max(partial_cor[,"value"],na.rm=TRUE);

for(i in 1:nrow(partial_cor) ){
	if(partial_cor[i,"Var1"]==partial_cor[i,"Var2"]){
		partial_cor[i,"value"]<-max_cor
	}
	
}

partial_cor[is.na(partial_cor[,"value"]),"value"]<-NA

########################################################################################################################

#clip_exon_binding_cor<-clip_binding_cor_spearman[clip_names_exon,clip_names_exon];
#clip_exon_binding_cor_melt<-melt(clip_exon_binding_cor);


## Edit body of pheatmap:::draw_colnames, customizing it to your liking
draw_colnames_45 <- function (coln, ...) {
    m = length(coln)
    x = (1:m)/m - 1/2/m
    grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
        hjust = 1, rot = 90, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
ns=asNamespace("pheatmap"))

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

#upper_tri_pearson<-get_upper_tri(clip_binding_cor_pearson);
#upper_tri_spearman<-get_upper_tri(clip_binding_cor_spearman);

melt_pearson<-melt(clip_binding_cor_pearson,na.rm=T);
melt_spearman<-melt(clip_binding_cor_spearman,na.rm=T);

write.table(melt_pearson, file="../result/RBP_exon_binding_cor_pearson.tsv",sep="\t",row.names=TRUE );
write.table(melt_spearman, file="../result/RBP_exon_binding_cor_spearman.tsv",sep="\t",row.names=TRUE );

source("multiplot.r")

#pheatmap(clip_binding_cor_pearson,cluster_rows=T,cluster_cols=T,fontsize_row=3,fontsize_col=3,main="Pearson correlation of RBP's binding sites");
#pheatmap(clip_binding_cor_spearman,cluster_rows=T,cluster_cols=T,fontsize_row=3,fontsize_col=3,main="Spearman correlation of RBP's binding sites");
p1<-ggplot(data = melt_pearson, aes(Var2, Var1, fill = value))+
 geom_tile()+
 scale_fill_gradient2(low = "white", high = "red", #mid = "yellow", 
   #midpoint = median(melt_pearson[,"value"]),
   limit = c(0,1), #space = "Lab", 
   name="Pearson\nCorrelation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 90, vjust = 1, 
    size = 0, hjust = 1),axis.text.y = element_text(angle = 0, vjust = 1, 
    size = 0, hjust = 1))+
 coord_fixed()+xlab("Pearson correlation")+ylab("RBP-region combinations")+guides(fill=FALSE)

p2<-ggplot(data = melt_spearman, aes(Var2, Var1, fill = value))+
 geom_tile()+
 scale_fill_gradient2(low = "white", high = "red", #mid = "yellow", 
   #midpoint = median(melt_spearman[,"value"]),
   limit = c(0,1),# space = "Lab", 
   name="Correlation\nCoefficient") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 90, vjust = 1, 
    size = 0, hjust = 1),axis.text.y = element_text(angle = 0, vjust = 1, 
    size = 0, hjust = 1))+
 coord_fixed()+
 xlab(paste0("<intron down region>   ","   <intron up region>    ","   <exon region>  "  ))+
 ggtitle("Spearman correlation")+ylab("RBP-region combinations")#+theme(legend.title=element_text("correlation"))#+guides(fill=FALSE)

p3<-ggplot(data = partial_cor, aes(Var2, Var1, fill = value))+
 geom_tile()+
 scale_fill_gradient2(low = "white", high = "red",# mid = "white", 
   #midpoint = median(partial_cor[,"value"]),
   limit = c(min(partial_cor[,"value"]),max(partial_cor[,"value"])),# space = "Lab", #min(partial_cor[,"value"]),max(partial_cor[,"value"])
   name="Correlation\nCoefficient") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 90, vjust = 1, 
    size = 0, hjust = 1),axis.text.y = element_text(angle = 0, vjust = 1, 
    size = 0, hjust = 1))+
 coord_fixed()+
 ggtitle("Partial Spearman correlation \nusing events' length as confounding factor")+
 xlab(paste0("<intron down region>   ","   <intron up region>    ","   <exon region>  "  ))+
 ylab("RBP-region combinations")#+theme(legend.title=element_text("correlation"))#+guides(fill=FALSE)

pdf("../result/Figure 1__RBP(exon_intron)_cor.pdf",height=6,width=13);
multiplot(p2,p3, cols=2,width=c(1,1) );
dev.off();


#clip_binding_cor_binary<-sapply(clip_binding_cor,function(x){if(x<0)return(-1);return(1); } );
#pdf("../anno/RBP(exon_intron)_cor_binary.pdf");
#pheatmap(clip_binding_cor_binary,cluster_rows=F,cluster_cols=F);
#dev.off();


#pdf("../result/RBP_attract_partial.pdf")

#dev.off()

