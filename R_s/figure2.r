library(ggplot2)
library(plyr)
library(reshape2)


plotFigure2<-function(){
  
  pdf("figure2_evaluate_each_feature.pdf",width=10,height=10);
  #jpeg("figure2.jpeg",width=800,height=800);
  
  plot.new();
  
  gl <- grid.layout(nrow=7, ncol=3, widths=c(1,1,1), heights=c(0.3,1,1,1,1,1,1) );
  
  vp.10 <- viewport(layout.pos.col=1, layout.pos.row=1);
  vp.20 <- viewport(layout.pos.col=2, layout.pos.row=1);
  vp.30 <- viewport(layout.pos.col=3, layout.pos.row=1);
    
  vp.11 <- viewport(layout.pos.col=1, layout.pos.row=1+1);
  vp.12 <- viewport(layout.pos.col=1, layout.pos.row=2+1);
  vp.13 <- viewport(layout.pos.col=1, layout.pos.row=3+1);
  vp.14 <- viewport(layout.pos.col=1, layout.pos.row=4+1);
  vp.15 <- viewport(layout.pos.col=1, layout.pos.row=5+1);
  vp.16 <- viewport(layout.pos.col=1, layout.pos.row=6+1);
  
  vp.21 <- viewport(layout.pos.col=2, layout.pos.row=1+1);
  vp.22 <- viewport(layout.pos.col=2, layout.pos.row=2+1);
  vp.23 <- viewport(layout.pos.col=2, layout.pos.row=3+1);
  vp.24 <- viewport(layout.pos.col=2, layout.pos.row=4+1);
  vp.25 <- viewport(layout.pos.col=2, layout.pos.row=5+1);
  vp.26 <- viewport(layout.pos.col=2, layout.pos.row=6+1);
  
  vp.31 <- viewport(layout.pos.col=3, layout.pos.row=1+1);
  vp.32 <- viewport(layout.pos.col=3, layout.pos.row=2+1);
  vp.33 <- viewport(layout.pos.col=3, layout.pos.row=3+1);
  vp.34 <- viewport(layout.pos.col=3, layout.pos.row=4+1);
  vp.35 <- viewport(layout.pos.col=3, layout.pos.row=5+1);
  vp.36 <- viewport(layout.pos.col=3, layout.pos.row=6+1);
  
  vps<-list();
  vps[[1]]<-vp.11;
  vps[[2]]<-vp.12;
  vps[[3]]<-vp.13;
  vps[[4]]<-vp.14;
  vps[[5]]<-vp.15;
  vps[[6]]<-vp.16;
  
  vps[[7]]<-vp.21;
  vps[[8]]<-vp.22;
  vps[[9]]<-vp.23;
  vps[[10]]<-vp.24;
  vps[[11]]<-vp.25;
  vps[[12]]<-vp.26;
  
  vps[[13]]<-vp.31;
  vps[[14]]<-vp.32;
  vps[[15]]<-vp.33;
  vps[[16]]<-vp.34;
  vps[[17]]<-vp.35;
  vps[[18]]<-vp.36;
  
  pushViewport(viewport(layout=gl)  );
  
  #pushViewport(vp.1);
  
  pushViewport(vp.10);  
  par(new=TRUE, fig=gridFIG(),mar=c(0.2,0.2,0.2,0.2));  
  
  textplot("A (density distribution)",cex=1);   
  #plot(c(1,2,3,4));
  popViewport();
  
  pushViewport(vp.20);
  par(new=TRUE, fig=gridFIG(),mar=c(0.2,0.2,0.2,0.2));  
  
  textplot("B (cpd)" ,cex=1); 
  #plot(c(1,2,3,4));
  popViewport();
  
  pushViewport(vp.30);
  par(new=TRUE, fig=gridFIG(),mar=c(0.2,0.2,0.2,0.2));  
  
  textplot("C (KS-test)",cex=1);  
  #plot(c(1,2,3,4));
  popViewport();
  
  goodFeatures<-alldatasetNA[,c("label","phylop","ss_3","asa_1","disorder_2","pfam","ptm") ];
  meltFeatures<-melt(goodFeatures,factorsAsStrings=FALSE);
  meltFeatures[,"variable"]<-as.character(meltFeatures[,"variable"])
  
  meltFeatures[,"variable"]<-factor( 
    abv(meltFeatures[,"variable"]), 
    levels=c(abv("phylop"),abv("ss_3"),abv("asa_1"),
             abv("disorder_2"),abv("pfam"),abv("ptm")
  ));
  
  cbbPalette <- c("#000000", "#AAAAAA","#000000","#E69F00", 
                  "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                  "#D55E00", "#CC79A7");
  
  for(i in 2:ncol(goodFeatures)){
    
    pushViewport(vps[[i-1]]);    
    par(new=TRUE, fig=gridFIG(),mar=c(2,2,2,1),cex=1);
    cat(paste(colnames(goodFeatures)[i],"\n" ,sep="") )
    densityNeutral<-density(goodFeatures[goodFeatures[,"label"]=="NEUTRAL",i],na.rm=TRUE );
    densityHGMD<-density(goodFeatures[goodFeatures[,"label"]=="HGMD",i],na.rm=TRUE );
    
    xMax<-max(c(densityNeutral$x,densityHGMD$x ) );xMin<-min(c(densityNeutral$x,densityHGMD$x ) );
    yMax<-max(c(densityNeutral$y,densityHGMD$y)  );yMin<-min(c(densityNeutral$y,densityHGMD$y)  );
    
    plot(densityNeutral$x,densityNeutral$y ,type="l",axes = F,
         ylim=c(yMin,yMax),xlim=c(xMin,xMax),col=cbbPalette[2]) ;
    
    lines(densityHGMD$x,densityHGMD$y,col=cbbPalette[1]) ;
    
    axis(1,cex.axis=0.7,mgp=c(0.5,0.5,0.5));
    axis(2,cex.axis=0.7,mgp=c(0.5,0.5,0.5));
    
    #plot(density(goodFeatures[goodFeatures[,"label"]=="HGMD",i],na.rm=TRUE ) ,add=TRUE) ;
    axisPos<-par("usr");
    start<-(axisPos[2]-axisPos[1])/20+axisPos[1];
    end<-axisPos[4]-(axisPos[4]-axisPos[3])/8;
    #text(start,end,labels=abv(colnames(goodFeatures)[i] ) ,cex=0.7,pos=4);
    mtext(text=abv(colnames(goodFeatures)[i] ),sid=3,cex=0.7,pos=4);
    
    popViewport();
    
  }
  
  phylopDataset<-alldatasetNA[,c("label","phylop")];
  ssDataset<-alldatasetNA[,c("label","ss_1","ss_2","ss_3","ss_4","ss_5","ss_6",
                           "ss_7","ss_8","ss_9","ss_10","ss_11","ss_12")];
  asaDataset<-alldatasetNA[,c("label","asa_1","asa_2","asa_3") ];
  disorderDataset<-alldatasetNA[,c("label","disorder_1","disorder_2",
                                 "disorder_3","disorder_4","disorder_5",
                                 "disorder_6","disorder_7","disorder_8",
                                 "disorder_9","disorder_10",
                                 "disorder_11","disorder_12")];
  
  pfamDataset<-alldatasetNA[,c("label","pfam")];
  ptmDataset<-alldatasetNA[,c("label","ptm")];
  
  
  features<-list();
  features[["phylop"]]<-phylopDataset;
  features[["ss"]]<-ssDataset;
  features[["asa"]]<-asaDataset;
  features[["disorder"]]<-disorderDataset;
  features[["pfam"]]<-pfamDataset;
  features[["ptm"]]<-ptmDataset;
  
  for(i in 2:ncol(goodFeatures) ){
    
    pushViewport(vps[[i-1+6]]);
    par(new=TRUE, fig=gridFIG(),mar=c(2,3,1,0) );
    #plot(c(1,2,3,4));
    
    plot(ecdf(goodFeatures[goodFeatures[,1]=="HGMD",i]  ) ,col=cbbPalette[1],
         pch=20,cex=0.5,axes = F,main="",ylim=c(0,1.4) );
    
    plot(ecdf(goodFeatures[goodFeatures[,1]=="NEUTRAL",i]  ) ,add=TRUE ,
         col=cbbPalette[2] ,pch=20,cex=0.5 ,axes = F,ylim=c(0,1.4) );
    
    axis(side=1, cex.axis=0.7, mgp=c(0.5,0.5,0.5) );
    axis(side=2, cex.axis=0.7, mgp=c(0.5,0.5,0.5) );
    
    axisPos<-par("usr");
    start<-(axisPos[2]-axisPos[1])/20+axisPos[1];
    text(start,axisPos[4]-0.25,labels=abv(names(goodFeatures)[i] ) ,cex=0.7,pos=4);
    #mtext(text=abv(names(goodFeatures)[i] ) ,cex=0.7,pos=4,side=3);
    
    popViewport();
    
  }
  
  #par( mcol=c(1,6) );
  #par(mar=c(0.5,0.5,0.5,0.5))
  for(i in 1:length(features) ){
    currentFeatures<-features[[i]];
    pValues<-c();
    
    for(j in 2:ncol(currentFeatures)){
      p<-ks.test(currentFeatures[currentFeatures[,1]=="NEUTRAL",j],
                 currentFeatures[currentFeatures[,1]=="HGMD",j])$p.value;
      
      if(p!=0){
        p<- (-1)*log10(p)
      }else{
        p<-16
        
      }
      p<-rnorm(1,p,0.4);
      pValues<-c(pValues,p);
      #pV
    }
    
    currentName<-names(features)[i];
    currentName<- abv(currentName);
    
    pushViewport(vps[[i+6+6]]);
    par(new=TRUE, fig=gridFIG(),mar=c(2,0,1,1));
    
    plot(x=pValues,y=rep(0.1,length(pValues)),ylim=c(-1,1)
         ,cex.main=0.5,axes = F,pch=4,xlim=c(-1,20),
         #ylab="",xlab="-log10( P value )",col="red",cex=0.7,mgp=c(2,1,.5));
         ylab="",xlab="",col=cbbPalette[3],cex=0.7,mgp=c(2,1,.5) );
    
    text(9,0.8,labels=(paste(currentName,"\n", " KS test P values -log10",sep="")),cex=0.7 );
    axis(side=1,pos=c(0,0.9),cex.axis=0.7 , mgp=c(1,0.35,0.35));
    
    if(i==length(features) ){
      par(xpd=TRUE);
      legend(8,-0.6,legend=c("hgmd","neutral"),col=c(cbbPalette[1],cbbPalette[2]),pch=c(15,15),border="white", );
      
    }
    
    popViewport();
    
  }
  #popViewport();
  #popViewport();
  
  #vpLabel <- viewport(layout.pos.col=2, layout.pos.row=7);
  #pushViewport(vpLabel);
  
  #par(new=TRUE, fig=gridFIG(),mar=c(1,1,1,1));
  #plot(c(0,1),c(0,1),yaxt="n",xaxt="n",xlab="",ylab="",type="l",axes=FALSE);
  #textplot("Neutral");
  #plot(c(0,1),c(0,1),yaxt="n",xaxt="n",xlab="",ylab="",type="l",col="#E69F00");
  #textplot("HGMD");
  #legend(0,0,legend=c("NEUTRAL","HGMD"),col=c("#000000","#E69F00") );
  #popViewport();
  
  dev.off();
  #ggplot(alldatasetNA[,c("label","ptm")] )+
  #  geom_histogram(aes(x=ptm,fill=label),binwidth=0.2,position="dodge" )+
  #  xlim(c(0,5));

}

setwd(workingDir)

plotFigure2();
