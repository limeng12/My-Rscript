library(reshape2);
library(ROCR);
#library(AUC)
library(gridBase);
library(grid);
library(pROC);

#setwd("/home/limeng/splicingSNP/R/"); 
#source("multiPlot.r"); 

plotFigure5<-function(){
  
  pdf("figure5_independent_test_and_Spanr_test.pdf",width=12,height=7);
  #jpeg("figure5.jpeg",width=1400,height=800);
  plot.new();
  
  # setup layout
  gl <- grid.layout(nrow=1, ncol=2);
  # grid.show.layout(gl)
  
  # setup viewports
  vp.1 <- viewport(layout.pos.col=1, layout.pos.row=1);
  vp.2 <- viewport(layout.pos.col=2, layout.pos.row=1);
  
  # init layout
  pushViewport(viewport(layout=gl));
  pushViewport(vp.1);
  
  par(new=TRUE, fig=gridFIG(),xpd=TRUE);
  #aucs<-sapply(featureAUC,"[",1);
  
  selectFeatures<-colnames(alldatasetTrain)[-1]; 
  alldatasetTrainSelect<-alldatasetTrain[,c("label",selectFeatures )];
  
  #for refseq
  modelrandomforest<-randomForest(formula=label~.,
                                  data=alldatasetTrainSelect,
                                  ntree=100,proximity=TRUE,
                                  replace=FALSE,nodesize=19,mtry=12);
  
  #for ensembl 
  #modelrandomforest<-randomForest(formula=label~., 
  #                                data=alldatasetTrainSelect, 
  #                                ntree=200,proximity=TRUE, 
  #                                replace=FALSE,nodesize=20,mtry=10); 
  
  #importantFeatures<<-importance(modelrandomforest); 
  #modelrandomforest<-ksvm(x=label~., data=alldatasetTrainSelect,prob.model = TRUE ); 
  #allDatasetTestROC<- 
  
  allDatasetTestSelect<-alldatasetTest[,c("label",selectFeatures )]; 
  predictValues<-predict(modelrandomforest,allDatasetTestSelect,type="prob")[,1];
  
  result<-data.frame(values=predictValues,label=allDatasetTestSelect$label );
  qplot(x=result$value,geom="density",data=result,color=result$label);
  
  cat(paste("true positive rate for 0.9 threshold:") );
  cat(sum(predictValues[allDatasetTestSelect$label=="HGMD"]>0.95)/sum(allDatasetTestSelect$label=="HGMD")  );
  
  tprFun<-function(x){
    t<-sum(predictValues[allDatasetTestSelect$label=="HGMD"]>x)/sum(allDatasetTestSelect$label=="HGMD");
    return(t);
  }
  
  fprFun<-function(x){
    t<-sum(predictValues[allDatasetTestSelect$label=="NEUTRAL"]>x)/sum(allDatasetTestSelect$label=="NEUTRAL");
    return(t);
  }
  
  f1Fun<-function(x){
    tp<-sum(predictValues[allDatasetTestSelect$label=="HGMD"]>x );   #tp
    fp<-sum(predictValues[allDatasetTestSelect$label=="NEUTRAL"]>x); #fp
    fn<-sum(predictValues[allDatasetTestSelect$label=="HGMD"]<x);    #fn
    f1<-(2*tp)/(2*tp+fp+fn);
    
    return(f1);
  }
  
  fdrFun<-function(x){
    t<-sum(predictValues[allDatasetTestSelect$label=="NEUTRAL"]>x)/sum(predictValues>x); #fdr
  }
  
  mccFun<-function(x){
    tp<-sum(predictValues[allDatasetTestSelect$label=="HGMD"]>x );      #tp
    tn<-sum(predictValues[allDatasetTestSelect$label=="NEUTRAL"]<x );   #tp
    fp<-sum(predictValues[allDatasetTestSelect$label=="NEUTRAL"]>x);    #fp
    fn<-sum(predictValues[allDatasetTestSelect$label=="HGMD"]<x);       #fn
    #print( paste0(tp," ",tn," ",fp," ",fn,"\n") );
    
    mcc<-( (tp*tn)-(fp*fn) )/( sqrt(tp+fp)*sqrt(tp+fn)*sqrt(tn+fp)*sqrt(tn+fn) ); 
    return(mcc);
  }
  
  thresholds<-seq(-0,1.1,0.0001);
  
  tpr<-sapply(thresholds,tprFun);
  fpr<-sapply(thresholds,fprFun);
  f1<- sapply(thresholds,f1Fun);
  fdr<-sapply(thresholds,fdrFun);
  mcc<-sapply(thresholds,mccFun);
  
  #troc<-roc(allDatasetTestSelect$label,predictValues );
  #plot(troc);
  
  #aucValue<-troc$auc;
  aucValue<-mean(sample(predictValues[allDatasetTestSelect$label=="HGMD"],1000000,replace=T) >
                   sample(predictValues[allDatasetTestSelect$label=="NEUTRAL"],1000000,replace=T));
    
  cutoffs <- data.frame(cut=thresholds, tpr=tpr, fpr=fpr,f1=f1,mcc=mcc);
  
  write.csv(cutoffs,file="cutoff.csv",quote=FALSE,row.names=FALSE);
  #plot(troc);
  
  jpeg("supplement_fig2.jpeg",width=900,height=600);
  cutoffs2<-melt(cutoffs,"cut");
  cuttpfp<-ggplot(cutoffs2) + geom_line(aes(x=cut,y=value,color=variable) )+
    geom_segment(aes(x = 0, y = 0.1, xend = 1, yend = 0.1) )+
    geom_segment(aes(x = 0.86, y = 0, xend = 0.86, yend = 1) )+
    geom_segment(aes(x = 0, y = 0.4765, xend = 1, yend = 0.4765) )+
    geom_text(aes(x=0 ,y=0.15, label="y=0.1") )+
    xlim(0,1)+ylim(0,1)+
    xlab("cutoff")+theme_classic();
  
  #print(cuttpfp);
  dev.off();
  
  #plot(troc,axes=FALSE);
  plot(x=fpr,y=tpr,xlab="False Positive Rate",ylab="True Positive Rate",type="l",xlim=c(0,1),ylim=c(0,1) );
  #abline(0,1,xlim=c(0,1),ylim=c(0,1) );
  segments(0,0,1,1);
  text(0.5,1.1,labels=paste("Test on independent test set\nAUC=", format(aucValue,digits=4),"\n",sep=" ") );
  
  #jpeg("tmp.jpeg")
  #plot(x=fpr,y=tpr,xlab="False Positive Rate",ylab="True Positive Rate",type="l",xlim=c(0,1),ylim=c(0,1) );
  #segments(0,0,1,1);
  #text(0.5,1.1,labels=paste("Test on independent test set\nAUC=", format(aucValue,digits=4),"\n",sep=" ") );
  
  #dev.off();
  #axis(side=1,at=seq(0,1,0.1),labels=seq(1,0,-0.1));
  #axis(side=2,at=seq(0,1,0.1),labels=seq(0,1,0.1));
  
  popViewport();
  pushViewport(vp.2);
  
  par(new=TRUE, fig=gridFIG());
  
  #psiPercent<-c(0.1309,0.19047,0.1978);
  psiPercent<-c(0.143,0.170,0.185);
  names(psiPercent)<-c("|PSI|>20%","10%<|PSI|<20%","|PSI|<10%");
  psiData<-data.frame(percent=psiPercent,names=names(psiPercent));
  
  #colnames(psiTestMelt)<-c("variable","value","percent")
  
  
  #b<-ggplot(psiData)+
  #  geom_bar(aes(x=names,y=percent),stat="identity",position="dodge",width=.7)+
  #  scale_x_discrete(limits=c("|PSI|>20","10<|PSI|<20","|PSI|<10"),
  #                      labels = c("|psi|>20", "20>|psi|>10","|psi|<10")  )+
  #  xlab("")+
  #  ylab("percentile")+
  #  ggtitle("Percentile of high disease score(>0.9) in different PSI group");
  #print(b, newpage = FALSE);
  barplot(psiPercent, 
          main="Percentile of high disease score(>0.9) in different PSI group",
          space=1.5,
          #horiz=TRUE,
          names.arg=names(psiPercent)  );
  
  popViewport();
  
  #multiplot(a,b,cols=2);
  
  dev.off();
  
}

setwd(workingDir)
plotFigure5();
