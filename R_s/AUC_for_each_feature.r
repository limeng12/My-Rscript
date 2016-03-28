
setwd(workingDir);

plotEachfeature<-function(){
  
  pdf("AUC.for.each.feature.pdf");
  
  text<-str_c("number of neutral exons:",nrow(alldataset[alldataset[,1]=="NEUTRAL",]),
              "\nnumber of HGMD exons:",nrow(alldataset[alldataset[,1]=="HGMD",])             
  );
  
  #textplot(text);
  
  
  #hist(normalData[,1],main="phylop score distribution",xlab="phylop score");
  #label<-c(rep("neutral",nrow(normalPhylop)),rep("hgmd",nrow(hgmdPhylop))  );
  #phylop<-rbind(normalPhylop[,1:2],hgmdPhylop[,1:2]);
  #phylop<-cbind(phylop,label);
  #phylop<-as.data.frame(phylop);
  #colnames(phylop)<-c("id","value","label");
  #phylop[,"label"]<-as.factor(phylop[,"label"]);
  
  #p<-ggplot(phylop)+geom_density(aes(x=value,color=label));
  #p<-p+ggtitle("phylop score between HGMD and neutral");
  #print(p);
  
  
  for(i in 2:ncol(alldataset)){
    
    
    featureName=colnames(alldataset)[i];
    currentFeature<-data.frame(alldataset[,1],alldataset[,i]);
    colnames(currentFeature)<-c("label","feature" );
    
    plot.new() ;
    
    # setup layout
    gl <- grid.layout(nrow=2, ncol=2);
    # grid.show.layout(gl)
    
    # setup viewports
    vp.1 <- viewport(layout.pos.col=1, layout.pos.row=1);
    vp.2 <- viewport(layout.pos.col=2, layout.pos.row=1);
    vp.3 <- viewport(layout.pos.col=c(1,2), layout.pos.row=2);
    
    # init layout
    pushViewport(viewport(layout=gl));
    
    pushViewport(vp.1);
    
    par(new=TRUE, fig=gridFIG());
    
    #plot(x = 1:10, y = 10:1)
    aroc<-1
    
    
    aroc<-roc(currentFeature[,1],currentFeature[,2],direction="auto");
    
    #if(21==i||30==i)
    #  aroc<-roc(currentFeature[,1],currentFeature[,2],direction=">")
    
    plot(aroc,main=paste("auc=",aroc$auc,sep="") );
    
    # done with the first viewport
    popViewport();
    
    pushViewport(vp.2);
    
    par(new=TRUE, fig=gridFIG());
    
    #plot(x = 1:10, y = 10:1)
    qqplot(currentFeature[currentFeature[,1]=="NEUTRAL",2],
           currentFeature[currentFeature[,1]=="HGMD",2]
           ,xlab="Neutral",ylab="HGMD",
           main=str_c(" KS test Pvalue",
           format(ks.test(currentFeature[currentFeature[,1]=="NEUTRAL",2],
                          currentFeature[currentFeature[,1]=="HGMD",2])$p.value,digits=2 )));
    
    abline(0,1);
    # done with the first viewport
    popViewport();
    
    # move to the next viewport
    pushViewport(vp.3)
    
    ggplotted <- qplot(x=1:10,y=10:1, 'point');
    # print our ggplot graphics here
    print(ggplotted, newpage = FALSE);
    
    # done with this viewport
    #popViewport(1)
    #currentFeature[,1]=factor(currentFeature[,1],labels=c("Neutral","HGMD"));
    
    
    a<-ggplot(currentFeature,mapping=aes(x=feature))+
      geom_density(mapping=aes(color=label),adjust=0.4)+
      geom_histogram(mapping=aes(fill=label,y=..density..),position="dodge")+
      xlab(featureName)+theme(legend.title=element_blank())+theme_classic()
    
    
    print(a,newpage = FALSE);
    #multiplot(a,b,c,cols=2)
    
    popViewport(2);
    
  }
  
  dev.off();

}

plotEachfeature();

