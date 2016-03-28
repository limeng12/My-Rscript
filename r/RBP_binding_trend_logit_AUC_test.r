library(ROCR)
library(ggplot2)
library(randomForest)
library(dplyr)
library(plyr)
#library("e1071")
setwd("E:\\limeng\\RBP_network\\r");
set.seed(100)
clip_binding_psi_data<-read.table(file="../anno/clip_binding_psi_var_median.tsv",sep="\t",header=TRUE,as.is=TRUE);
#fileter_psi_names<-colnames(RBP_binding_psi_data)[ grepl("EN",colnames(RBP_binding_psi_data) )];
#clip_binding_psi_data<-clip_binding_psi_data[,setdiff(colnames(clip_binding_psi_data),c("start.y","end.y","start","end",fileter_psi_names) ) ];

#filter
clip_binding_psi_data_filter<-subset(clip_binding_psi_data, (median_psi>0.9|median_psi<0.1)&var_psi<0.1&count_psi>70 );

is_constitute<-sapply(clip_binding_psi_data_filter[,"median_psi"],function(x){if(x>0.9) return(TRUE);return(FALSE) } );

clip_binding_psi_data_filter_class<-cbind(is_constitute,clip_binding_psi_data_filter);

clip_binding_psi_data_filter_class<-clip_binding_psi_data_filter_class


clip_binding_psi_data_filter_class_as<-clip_binding_psi_data_filter_class[clip_binding_psi_data_filter_class[,1]==FALSE,]
clip_binding_psi_data_filter_class_constitute<-sample_n(clip_binding_psi_data_filter_class[clip_binding_psi_data_filter_class[,1]==TRUE,],nrow(clip_binding_psi_data_filter_class_as))

clip_binding_psi_data_filter_class_equal_redun<-rbind( clip_binding_psi_data_filter_class_constitute,clip_binding_psi_data_filter_class_as);
write.table(clip_binding_psi_data_filter_class_equal_redun,file="../anno/RBP_binding_psi_for_logistic_regression.tsv",sep="\t",quote=FALSE,row.names=FALSE);

################################################################################overlap with UCSC knownAlt hg19############################################################
clip_binding_psi_data_filter_class_equal<-clip_binding_psi_data_filter_class_equal_redun[,-c(2:6)];

#ucsc_known_alt<-read.table("../anno/ucsc_hg19_knownAlt.bed",header=FALSE,as.is=TRUE,sep=",");
#colnames(ucsc_known_alt)<-c("chr","start","end","name","score","strand");

ensembl_known_alt<-read.table("../anno/ensembl_alt.csv",header=TRUE,as.is=TRUE,sep=",");
colnames(ensembl_known_alt)<-c("chr","start","end");
ensembl_known_alt[,"chr"]<-paste0("chr",ensembl_known_alt[,"chr"]);
ensembl_known_alt_unique<-unique(ensembl_known_alt);

center_exons_split<-strsplit(clip_binding_psi_data_filter_class_equal_redun[,"exonId"],":")

center_exons_start<-as.numeric(sapply(center_exons_split,"[",5));
center_exons_end<-as.numeric(sapply(center_exons_split,"[",6));
center_exons_chr<-as.character(sapply(center_exons_split,"[",1))


clip_binding_psi_data_filter_class_equal_redun<-cbind(clip_binding_psi_data_filter_class_equal_redun,center_exons_start,center_exons_end,center_exons_chr);
number_of_unique_center_exons<-sum(!duplicated(clip_binding_psi_data_filter_class_equal_redun[,c("center_exons_start","center_exons_end","center_exons_chr")]));
print(paste0("Number of unique center exons contain by the events: ",number_of_unique_center_exons))

clip_binding_psi_data_filter_class_equal_redun[,"center_exons_start"]<-as.numeric(clip_binding_psi_data_filter_class_equal_redun[,"center_exons_start"]);
clip_binding_psi_data_filter_class_equal_redun[,"center_exons_end"]<-as.numeric(clip_binding_psi_data_filter_class_equal_redun[,"center_exons_end"]);
clip_binding_psi_data_filter_class_equal_redun[,"center_exons_chr"]<-as.character(clip_binding_psi_data_filter_class_equal_redun[,"center_exons_chr"]);


clip_binding_psi_data_filter_class_equal_redun_knownAlt_join<-inner_join(clip_binding_psi_data_filter_class_equal_redun,ensembl_known_alt,by=c("center_exons_chr"="chr","center_exons_start"="start","center_exons_end"="end"))

clip_binding_psi_data_filter_class_equal_redun_knownAlt_join_unique<-clip_binding_psi_data_filter_class_equal_redun_knownAlt_join[!duplicated(clip_binding_psi_data_filter_class_equal_redun_knownAlt_join[,"exonId"]),]

######################################################################################################################################################################################


#####################################################################################PSI confidence interval estimation###############################################################


######################################################################################################################################################################################

clip_binding_psi_data_filter_class<-clip_binding_psi_data_filter_class_equal;

#set.seed(12132)
clip_binding_psi_data_filter_class<-clip_binding_psi_data_filter_class[sample(1:nrow(clip_binding_psi_data_filter_class),nrow(clip_binding_psi_data_filter_class) ),];

write.table(clip_binding_psi_data_filter_class,file="../anno/clip_binding_psi_data_class.tsv",sep="\t",quote=FALSE,row.names=FALSE);

print(paste0("number of class: ",nrow(clip_binding_psi_data_filter_class) ));
line3_2<-as.integer(nrow(clip_binding_psi_data_filter_class)/3*2);

#clip_binding_psi_data_filter_class_logitModel<-glm(is_constitute~.,unique(clip_binding_psi_data_filter_class[1:line3_2,]), family = "binomial",control = list(maxit = 5000),trace = TRUE, epsilon = 1e-14);
#clip_binding_psi_data_filter_class_logitModel<-naiveBayes(is_constitute~.,clip_binding_psi_data_filter_class[1:line3_2,] );

clip_binding_psi_data_filter_class[,"is_constitute"]<-as.factor(clip_binding_psi_data_filter_class[,"is_constitute"]);
clip_binding_psi_data_filter_class_logitModel<-randomForest(is_constitute~.,clip_binding_psi_data_filter_class[1:line3_2,],ntree=100,mtry=12);

clip_binding_psi_data_filter_class_test<-clip_binding_psi_data_filter_class[(line3_2+1):nrow(clip_binding_psi_data_filter_class),];

#clip_binding_psi_data_filter_class_test_pre<-predict(clip_binding_psi_data_filter_class_logitModel,clip_binding_psi_data_filter_class_test,type='response');
clip_binding_psi_data_filter_class_test_pre<-predict(clip_binding_psi_data_filter_class_logitModel,clip_binding_psi_data_filter_class_test,type='prob')[,2];
#clip_binding_psi_data_filter_class_test_pre<-predict(clip_binding_psi_data_filter_class_logitModel,clip_binding_psi_data_filter_class_test,type='raw')[,2];


predict_psi<-cbind(clip_binding_psi_data_filter_class_test_pre,clip_binding_psi_data_filter_class_test);
#write.table(predict_psi,file="../result/logistic regression predict.tsv",sep="\t",row.names=FALSE)


clip_binding_psi_data_filter_class_test_pred<-prediction(clip_binding_psi_data_filter_class_test_pre,clip_binding_psi_data_filter_class_test[,"is_constitute"]);
prf_counts<-performance(clip_binding_psi_data_filter_class_test_pred,measure="tpr",x.measure="fpr");
auc<-performance(clip_binding_psi_data_filter_class_test_pred,measure="auc");
aucv_counts<-signif(auc@y.values[[1]] ,digits=3);

#clip_binding_psi_data_filter_class_logitModel_sum<-summary(clip_binding_psi_data_filter_class_logitModel);
#coe_logit<-clip_binding_psi_data_filter_class_logitModel_sum$coefficients[order(clip_binding_psi_data_filter_class_logitModel_sum$coefficients[,4]),]

#logit_coeffcients<-coe_logit[1:20,1]
#logit_coeffcients<-signif(logit_coeffcients,digits=2)
#logit_p_value<-coe_logit[1:20,4]
#logit_p_value<--log10(logit_p_value+ (1e-16) );

#pdf("../anno/Logistic regression Coeffcients and p-value.pdf",width=10,height=7);


#logit_p_value_format_names<-names(logit_p_value_format);
#logit_p_value_format_names<-factor(logit_p_value_format_names,levels=logit_p_value_format_names );

#logit_p_value_format<-signif(logit_p_value,digits=4);
#p<-ggplot()+geom_bar(aes(x=logit_p_value_format_names,y=logit_p_value_format),stat="identity")+
# geom_text(aes(x=logit_p_value_format_names,y=logit_p_value_format,label=logit_coeffcients),vjust=-0.2,size=3)+theme_classic()+
# theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+ylab("p-value(-log10)")+ggtitle("top 20 P-values and coefficients for Logistic regression")+
# xlab("RBP");
#print(p);

#dev.off();



pdf("../result/AUC of logistic model counts.pdf");
plot(prf_counts,main=paste0("Logistic regression of RBP and PSI, AUC = ",aucv_counts) );
dev.off();

fitted_results <- ifelse(clip_binding_psi_data_filter_class_test_pre > 0.5,TRUE,FALSE);

clip_binding_psi_data_filter_class_test_result<-cbind(fitted_results,clip_binding_psi_data_filter_class_test);

write.table(clip_binding_psi_data_filter_class_test_result,file="../result/Logistic_regression_of_RBP_counts.tsv",sep="\t");
ACU<-(1-with(clip_binding_psi_data_filter_class_test_result,mean(fitted_results!=is_constitute)) );

###########################################################################################################################################################
#Binary the number of RBP in each site

clip_binding_psi_data_filter_class_binary<-clip_binding_psi_data_filter_class;

clip_binding_psi_data_filter_class_binary[,2:ncol(clip_binding_psi_data_filter_class_binary)]<-apply(clip_binding_psi_data_filter_class[,-1],c(1,2),function(x){if(x>0)return(1);return(0) } );
clip_binding_psi_data_filter_class_binary<-as.data.frame(clip_binding_psi_data_filter_class_binary)
#clip_binding_psi_data_filter_class_binary_logitModel<-glm(is_constitute~.,clip_binding_psi_data_filter_class_binary[1:line3_2,], family = "binomial",control = list(maxit = 500));

clip_binding_psi_data_filter_class_binary[,"is_constitute"]<-as.factor(clip_binding_psi_data_filter_class_binary[,"is_constitute"]);
clip_binding_psi_data_filter_class_binary_logitModel<-randomForest(is_constitute~.,clip_binding_psi_data_filter_class_binary[1:line3_2,],ntree=100,mtry=12);

clip_binding_psi_data_filter_class_binary_test<-clip_binding_psi_data_filter_class_binary[(line3_2+1):nrow(clip_binding_psi_data_filter_class_binary),];

#clip_binding_psi_data_filter_class_binary_test_pre<-predict(clip_binding_psi_data_filter_class_binary_logitModel,clip_binding_psi_data_filter_class_binary_test,type='response');
clip_binding_psi_data_filter_class_binary_test_pre<-predict(clip_binding_psi_data_filter_class_binary_logitModel,clip_binding_psi_data_filter_class_binary_test,type='prob')[,2];

clip_binding_psi_data_filter_class_binary_test_pred<-prediction(clip_binding_psi_data_filter_class_binary_test_pre,clip_binding_psi_data_filter_class_binary_test[,"is_constitute"]);
    
prf_binary<-performance(clip_binding_psi_data_filter_class_binary_test_pred,measure="tpr",x.measure="fpr");
auc<-performance(clip_binding_psi_data_filter_class_binary_test_pred,measure="auc");
aucv_binary<-signif(auc@y.values[[1]] ,digits=3);

pdf("../result/AUC of logistic model binary.pdf");
plot(prf_binary,main=paste0("Logistic regression of RBP binary and PSI, AUC = ",aucv_binary) );
dev.off();

fitted_results <- ifelse(clip_binding_psi_data_filter_class_binary_test_pre > 0.5,TRUE,FALSE);
clip_binding_psi_data_filter_class_binary_test_result<-cbind(fitted_results,clip_binding_psi_data_filter_class_binary_test);
write.table(clip_binding_psi_data_filter_class_binary_test_result,file="../result/Logistic_regression_of_RBP_binary.tsv",sep="\t");


################################################################################################################################################################
#use RBP's expression level

gene_expression<-read.table("../anno/gene_expression_table.tsv",header=TRUE,as.is=TRUE,);

median_exp<-apply(gene_expression[,-1],1,median);
variance_exp<-apply(gene_expression[,-1],1,var);

gene_expression_sum<-data.frame(gene=gene_expression[,1],median_exp=median_exp,variance_exp=variance_exp,stringsAsFactors=FALSE);

clip_binding_psi_data_filter_class_exp<-clip_binding_psi_data_filter_class_binary

RBP_names<-sapply(strsplit(colnames(clip_binding_psi_data_filter_class_exp),"_" ),"[",1 );
for(i in 2:ncol(clip_binding_psi_data_filter_class_exp) ){
  
  rbp_exp<-rep(gene_expression_sum[gene_expression_sum[,1]==RBP_names[i],2],nrow(clip_binding_psi_data_filter_class_exp) );
  rbp_exp<-rbp_exp*clip_binding_psi_data_filter_class_exp[,i];
  
  clip_binding_psi_data_filter_class_exp[,i]<-rbp_exp
  
}

gene_expression_sum_rbp<-gene_expression_sum[is.element(gene_expression_sum[,1],RBP_names),]; 
write.table(gene_expression_sum_rbp,file="../result/RBP expresion variance.tsv",sep="\t",quote=FALSE)


#gene_expression_sum_rbp_var_small<-gene_expression_sum_rbp[gene_expression_sum_rbp[,3]<100,];

#gene_expression_sum_rbp_var_small_exon<-paste0(gene_expression_sum_rbp_var_small[,1],"_exon")
#gene_expression_sum_rbp_var_small_intron_down<-paste0(gene_expression_sum_rbp_var_small[,1],"_intron_down")
#gene_expression_sum_rbp_var_small_intron_up<-paste0(gene_expression_sum_rbp_var_small[,1],"_intron_up")

#gene_expression_sum_rbp_var_small_names<-c(gene_expression_sum_rbp_var_small_exon,gene_expression_sum_rbp_var_small_intron_down,gene_expression_sum_rbp_var_small_intron_up )

#clip_binding_psi_data_filter_class_exp<-clip_binding_psi_data_filter_class_exp[,c("is_constitute",gene_expression_sum_rbp_var_small_names ) ];

clip_binding_psi_data_filter_class_exp[,"is_constitute"]<-as.factor(clip_binding_psi_data_filter_class_exp[,"is_constitute"]);
#clip_binding_psi_data_filter_class_exp_logitModel<-glm(is_constitute~.,clip_binding_psi_data_filter_class_exp[1:line3_2,], family = "binomial",control = list(maxit = 500));
clip_binding_psi_data_filter_class_exp_logitModel<-randomForest(is_constitute~.,clip_binding_psi_data_filter_class_exp[1:line3_2,]);


clip_binding_psi_data_filter_class_exp_test<-clip_binding_psi_data_filter_class_exp[(line3_2+1):nrow(clip_binding_psi_data_filter_class),];

#clip_binding_psi_data_filter_class_exp_test_pre<-predict(clip_binding_psi_data_filter_class_exp_logitModel,clip_binding_psi_data_filter_class_exp_test,type='response');
clip_binding_psi_data_filter_class_exp_test_pre<-predict(clip_binding_psi_data_filter_class_exp_logitModel,clip_binding_psi_data_filter_class_exp_test,type='prob')[,2];


clip_binding_psi_data_filter_class_exp_test_pred<-prediction(clip_binding_psi_data_filter_class_exp_test_pre,clip_binding_psi_data_filter_class_exp_test[,"is_constitute"]);
prf_expression<-performance(clip_binding_psi_data_filter_class_exp_test_pred,measure="tpr",x.measure="fpr");
auc<-performance(clip_binding_psi_data_filter_class_exp_test_pred,measure="auc");
aucv_exp<-signif(auc@y.values[[1]] ,digits=3);

pdf("../result/AUC of logistic model expression.pdf");
plot(prf_expression,main=paste0("Logistic regression of RBP expression and PSI, AUC = ",aucv_exp) );
dev.off();


fitted_results <- ifelse(clip_binding_psi_data_filter_class_exp_test_pre > 0.5,TRUE,FALSE);
clip_binding_psi_data_filter_class_exp_test_result<-cbind(fitted_results,clip_binding_psi_data_filter_class_exp_test);
write.table(clip_binding_psi_data_filter_class_exp_test_result,file="../result/Logistic_regression_of_RBP_exp.tsv",sep="\t");

randomForest_RBP_psi<-randomForest(is_constitute~.,clip_binding_psi_data_filter_class_exp,ntree=100,mtry=12);

impFeatures<-importance(randomForest_RBP_psi);

pdf("../result/feature importance.pdf",width=7,height=5);

#plot(impFeatures)
impFeatures_sort<-impFeatures[order(impFeatures,decreasing=TRUE),];
impFeatures_sort_names<-factor(names(impFeatures_sort),levels=names(impFeatures_sort))

barplot(impFeatures_sort[1:10],xlab="RBP names",ylab="importance by RandomForest");
#p_feature_importance<-ggplot()+geom_bar(aes(x=impFeatures_sort_names,y=impFeatures_sort ),stat="identity")+
  #geom_text(aes(x=names(impFeatures_sort),y=impFeatures_sort,label=logit_coeffcients),vjust=-0.2,size=3)
  #theme_classic()+
  #theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=3))+ylab("importance")+ggtitle("Feature importance by randomforest")+
  #xlab("RBP");
#print(p);

dev.off();



crossValidateRandomForest<-function(data,crossNum=10){
  allSampleIDs<-1:nrow(data);
  onePortionSize<-floor( length(allSampleIDs)/crossNum ); 
  
  sumRoc<-0;
  
  all_predict<-c();
  all_is_constitute<-c();
  
  for(i in 1:crossNum){
    set.seed(1000);
    testID<-sample(allSampleIDs,onePortionSize,replace=FALSE); 
    allSampleIDs<-setdiff(allSampleIDs,testID); 
    trainID<-setdiff(  1:nrow(data),  testID  ); 
    

    #model<-glm(is_constitute~.,data[trainID,], family = "binomial",control = list(maxit = 50));
    randomForest_model<-randomForest(is_constitute~.,data[trainID,],ntree=100,mtry=12);

    #troc<-(roc( (data[testID,])$label, predict(modelrandomforest,data[testID,],type="prob")[,1] )$auc) [1];
    predict<-predict(randomForest_model,data[testID,],type='prob')[,2]; 
    is_constitute<-(data[testID,])$is_constitute; 
    
    all_predict<-c(all_predict,predict) ; 
    all_is_constitute<-c(all_is_constitute,is_constitute ) ; 
  }
  
  #return(data.frame(prob=allProb,label=allLabel  )  ); 
  
  clip_binding_psi_data_filter_class_test_pred<-prediction(all_predict,all_is_constitute );
  #prf_counts<-performance(clip_binding_psi_data_filter_class_test_pred,measure="tpr",x.measure="fpr");
  return(clip_binding_psi_data_filter_class_test_pred)
}
#aucv_exp<-signif(auc@y.values[[1]] ,digits=3);

rbp_counts_pred_cv<-crossValidateRandomForest(clip_binding_psi_data_filter_class);
rbp_counts_pred_prf<-performance(rbp_counts_pred_cv,measure="tpr",x.measure="fpr");
rbp_counts_pred_auc<-signif(performance(rbp_counts_pred_cv,measure="auc")@y.values[[1]],digits=3);

rbp_binary_pred_cv<-crossValidateRandomForest(clip_binding_psi_data_filter_class_binary);
rbp_binary_pred_prf<-performance(rbp_binary_pred_cv,measure="tpr",x.measure="fpr");
rbp_binary_pred_auc<-signif(performance(rbp_binary_pred_cv,measure="auc")@y.values[[1]],digits=3);

rbp_expression_pred_cv<-crossValidateRandomForest(clip_binding_psi_data_filter_class_exp);
rbp_expression_pred_prf<-performance(rbp_expression_pred_cv,measure="tpr",x.measure="fpr");
rbp_expression_pred_auc<-signif(performance(rbp_expression_pred_cv,measure="auc")@y.values[[1]],digits=3);


pdf("../result/Figure 4 logistic region of RBP-position combinations with PSI.pdf")
par(mfrow=c(2,2));
par(mar=c(5.1,4.1,2.1,2.1))

plot(rbp_counts_pred_prf,main=paste0("(A)Classification of PSI\n using RBP counts, AUC = ",rbp_counts_pred_auc) );

plot(rbp_binary_pred_prf,main=paste0("(B)Classification of PSI\n using RBP binary, AUC = ",rbp_binary_pred_auc) );

plot(rbp_expression_pred_prf,main=paste0("(C)Classification of PSI\n using RBP expression, AUC = ",rbp_expression_pred_auc) );

par(mar=c(7.1,4.1,2.1,2.1))
barplot(impFeatures_sort[1:10],xlab="",main="(D)Top ten importance RBPs\n by RandomForest",cex.names=0.7,las=2);

dev.off();

