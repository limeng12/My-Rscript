#predictClinvar
#clinsig5PredictByPsi
setwd(RPath);
source("multiPlot.r");

plotFigure6<-function(){
  setwd(workingDir)
  pdf("figure6.pdf")
  
  comparePsi<-data.frame(
    value=c(predictClinvar5,clinsig5PredictByPsi),
    label=c(rep("model",length(predictClinvar5) ),rep("PSI",length(clinsig5PredictByPsi) )  )
    );
  
  p1<-ggplot( subset(comparePsi,label=="model" )  )+geom_histogram(aes(x=value))+
    xlab("disease cansing probability")+ggtitle("predict by our model(Pathogenic snp on splicing site)");
  
  p2<-ggplot( subset(comparePsi,label=="PSI" )  )+geom_histogram(aes(x=value))+
    xlab("|dPSI|")+ggtitle("predict by SPANR(Pathogenic snp on splicing site)");
  
  multiplot(p1,p2,cols=1);
  dev.off();
  
}


plotFigure6();


plotFigure7<-function(){
  setwd(workingDir);
  pdf("figure7.pdf");
  
  comparePsi<-data.frame(
    value=c(predictClinvar5,predictClinvar2),
    label=c(rep("our model on disease",length(predictClinvar5) ),rep("our model on benign",length(predictClinvar2) )  )
  );
  
  p1<-ggplot( subset(comparePsi,label=="our model on disease" )  )+geom_histogram(aes(x=value))+
    xlab("disease causing probability")+ggtitle("predict by our model(Pathogenic snp on splicing site)" );
  
  p2<-ggplot( subset(comparePsi,label=="our model on benign" )  )+geom_histogram(aes(x=value))+
    xlab("disease causing probability")+ggtitle("predict by our model(benign snp on splicing site)" );
  
  multiplot(p1,p2,cols=1);
  dev.off();
  
}

plotFigure7();

plotFigure8<-function(){
  setwd(workingDir);
  pdf("figure8.pdf");
  
  comparePsi<-data.frame(
    value=c(clinsig5PredictByPsi,clinsig2PredictByPsi),
    label=c(rep("SPANR on disease",length(clinsig5PredictByPsi) ),rep("SPANR on benign",length(clinsig2PredictByPsi) ) )  
  );
  
  p1<-ggplot( subset(comparePsi,label=="SPANR on disease" )  )+geom_histogram(aes(x=value))+
    xlab("|PSI|")+ggtitle("predict by SPANR(Pathogenic snp on splicing site)")
  
  p2<-ggplot( subset(comparePsi,label=="SPANR on benign" )  )+geom_histogram(aes(x=value))+
    xlab("|PSI|")+ggtitle("predict by SPANR(benign snp on splicing site)");
  
  multiplot(p1,p2,cols=1);
  dev.off();
  
}

plotFigure8();




