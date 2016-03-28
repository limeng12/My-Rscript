#alldataset
pdf("figure2.pdf",width=8,height=10)

plot.new() 
gl <- grid.layout(nrow=6, ncol=2);
pushViewport(viewport(layout=gl));

##################################phylop socre###############################################################
currentFeature=data.frame(label=alldataset[,1],value=alldataset[,2] ,stringsAsFactors =FALSE);
#format(ks.test(currentFeature[currentFeature[,1]==0,2],currentFeature[currentFeature[,1]==1,2])$p.value,digits=2 )
##phylop score histogram plot part
vp.1 <- viewport(layout.pos.col=1, layout.pos.row=1);
pushViewport(vp.1);

a<-ggplot(currentFeature,mapping=aes(x=currentFeature[,2]))+
  geom_density(mapping=aes(color=label),adjust=0.4)+
  xlab("phylop score")+theme(legend.title=element_blank());
print(a,newpage = FALSE);

popViewport();

#phylop score p-value plot part
vp.2 <- viewport(layout.pos.col=2, layout.pos.row=1);
pushViewport(vp.2);
#par(new=TRUE, fig=gridFIG(),mar=c(2.0,1.0,2.0,1.0))

p.value=format(ks.test(currentFeature[currentFeature[,1]=="NEUTRAL",2],
                       currentFeature[currentFeature[,1]=="HGMD",2])$p.value,digits=2 );
#plot(p.value,xlab="feature index",ylab="p value",main="Kolmogorov–Smirnov test p value")
b<-ggplot()+geom_point(aes(x=seq(1:length(p.value)),y=p.value))+
  xlab("feature index")+scale_x_continuous(breaks=seq(1:3))+
  ggtitle("Kolmogorov Smirnov test")+ylab("p value");

print(b,newpage = FALSE);
#plot(0.5,0.5);
popViewport();

###################################Secondary structure plot##################################################################
currentFeature=data.frame(label=alldataset[,1],value=alldataset[,5] ,stringsAsFactors =FALSE);
currentFeatures=data.frame(label=alldataset[,1],value=alldataset[,3:14] ,stringsAsFactors =FALSE)

vp.3 <- viewport(layout.pos.col=1, layout.pos.row=2);
pushViewport(vp.3);

a<-ggplot(currentFeature,mapping=aes(x=currentFeature[,2]))+
  geom_density(mapping=aes(color=label),adjust=0.4)+
  xlab("Secondary structure")+
  theme(legend.title=element_blank());
print(a,newpage = FALSE);

popViewport();


vp.4 <- viewport(layout.pos.col=2, layout.pos.row=2);
pushViewport(vp.4);
par(new=TRUE, fig=gridFIG(),mar=c(1.0,1.0,1.0,1.0))
p.values<-c();

for(i in 2:ncol(currentFeatures)){
  
  p.value=format(ks.test(currentFeatures[currentFeatures[,1]=="NEUTRAL",i],
                         currentFeatures[currentFeatures[,1]=="HGMD",i])$p.value,digits=2 );
  p.values[i-1]<-p.value;
  
}
#plot(p.values,xlab="feature index",ylab="p value",main="Kolmogorov–Smirnov test p value")
b<-ggplot()+geom_point(aes(x=seq(1:length(p.values)),y=p.values))+
  xlab("feature index")+scale_x_continuous(breaks=seq(1:12))+
  ggtitle("Kolmogorov Smirnov test")+ylab("p value");
print(b,newpage = FALSE);
popViewport();

###################################ASA plot##################################################################
currentFeature=data.frame(label=alldataset[,1],value=alldataset[,15] ,stringsAsFactors =FALSE);
currentFeatures=data.frame(label=alldataset[,1],value=alldataset[,15:17] ,stringsAsFactors =FALSE)

vp.5 <- viewport(layout.pos.col=1, layout.pos.row=3);
pushViewport(vp.5);

a<-ggplot(currentFeature,mapping=aes(x=currentFeature[,2]))+
  geom_density(mapping=aes(color=label),adjust=0.4)+
  xlab("ASA area")+theme(legend.title=element_blank());
print(a,newpage = FALSE);

popViewport();


vp.6 <- viewport(layout.pos.col=2, layout.pos.row=3);
pushViewport(vp.6);
par(new=TRUE, fig=gridFIG(),mar=c(1.0,1.0,1.0,1.0))
p.values<-c();

for(i in 2:ncol(currentFeatures)){
  
  p.value=format(ks.test(currentFeatures[currentFeatures[,1]=="NEUTRAL",i],currentFeatures[currentFeatures[,1]=="HGMD",i])$p.value,digits=2 );
  p.values[i-1]<-p.value;
  
}
#plot(p.values,xlab="feature index",ylab="p value",main="Kolmogorov–Smirnov test p value")
b<-ggplot()+geom_point(aes(x=seq(1:length(p.values)),y=p.values))+xlab("feature index")+scale_x_continuous(breaks=seq(1:12))+ggtitle("Kolmogorov Smirnov test")+ylab("p value");
print(b,newpage = FALSE);

popViewport();


###################################disorder plot##################################################################
currentFeature=data.frame(label=alldataset[,1],value=alldataset[,29] ,stringsAsFactors =FALSE);
currentFeatures=data.frame(label=alldataset[,1],value=alldataset[,18:29] ,stringsAsFactors =FALSE)

vp.7 <- viewport(layout.pos.col=1, layout.pos.row=4);
pushViewport(vp.7);

a<-ggplot(currentFeature,mapping=aes(x=currentFeature[,2]))+geom_density(mapping=aes(color=label),adjust=0.4)+
  xlab("disorder probability")+theme(legend.title=element_blank());
print(a,newpage = FALSE);

popViewport();


vp.8 <- viewport(layout.pos.col=2, layout.pos.row=4);
pushViewport(vp.8);
par(new=TRUE, fig=gridFIG(),mar=c(1.0,1.0,1.0,1.0))
p.values<-c();

for(i in 2:ncol(currentFeatures)){
  
  p.value=format(ks.test(currentFeatures[currentFeatures[,1]=="NEUTRAL",i],currentFeatures[currentFeatures[,1]=="HGMD",i])$p.value,digits=2 );
  p.values[i-1]<-p.value;
  
}
#plot(p.values,xlab="feature index",ylab="p value",main="Kolmogorov–Smirnov test p value");
b<-ggplot()+geom_point(aes(x=seq(1:length(p.values)),y=p.values))+xlab("feature index")+scale_x_continuous(breaks=seq(1:12))+ggtitle("Kolmogorov Smirnov test")+ylab("p value");
print(b,newpage = FALSE);
popViewport();


###################################pfam plot##################################################################
currentFeature=data.frame(label=alldataset[,1],value=alldataset[,30] ,stringsAsFactors =FALSE);
currentFeatures=data.frame(label=alldataset[,1],value=alldataset[,30] ,stringsAsFactors =FALSE)

vp.9 <- viewport(layout.pos.col=1, layout.pos.row=5);
pushViewport(vp.9);

a<-ggplot(currentFeature,mapping=aes(x=currentFeature[,2]))+geom_density(mapping=aes(color=label),adjust=0.4)+
  xlab("pfam coverage percentage")+theme(legend.title=element_blank());
print(a,newpage = FALSE);

popViewport();


vp.10 <- viewport(layout.pos.col=2, layout.pos.row=5);
pushViewport(vp.10);
par(new=TRUE, fig=gridFIG(),mar=c(1.0,1.0,1.0,1.0))
p.values<-c();

for(i in 2:ncol(currentFeatures)){
  
  p.value=format(ks.test(currentFeatures[currentFeatures[,1]=="NEUTRAL",i],currentFeatures[currentFeatures[,1]=="HGMD",i])$p.value,digits=2 );
  p.values[i-1]<-p.value;
  
}
#plot(p.values,xlab="feature index",ylab="p value",main="Kolmogorov–Smirnov test p value");
b<-ggplot()+geom_point(aes(x=seq(1:length(p.values)),y=p.values))+xlab("feature index")+scale_x_continuous(breaks=seq(1:3))+ggtitle("Kolmogorov Smirnov test")+ylab("p value");
print(b,newpage = FALSE);
popViewport();


###################################ptm plot##################################################################
currentFeature=data.frame(label=alldataset[,1],value=alldataset[,31] ,stringsAsFactors =FALSE);
currentFeatures=data.frame(label=alldataset[,1],value=alldataset[,31] ,stringsAsFactors =FALSE)

vp.11 <- viewport(layout.pos.col=1, layout.pos.row=6);
pushViewport(vp.11);

a<-ggplot(currentFeature,mapping=aes(x=currentFeature[,2]))+geom_density(mapping=aes(color=label),adjust=0.4)+
  xlab("ptm probability")+theme(legend.title=element_blank());
print(a,newpage = FALSE);

popViewport();


vp.12 <- viewport(layout.pos.col=2, layout.pos.row=6);
pushViewport(vp.12);
par(new=TRUE, fig=gridFIG(),mar=c(1.0,1.0,1.0,1.0))
p.values<-c();

for(i in 2:ncol(currentFeatures)){
  
  p.value=format(ks.test(currentFeatures[currentFeatures[,1]=="NEUTRAL",i],currentFeatures[currentFeatures[,1]=="HGMD",i])$p.value,digits=2 );
  p.values[i-1]<-p.value;
  
}
#plot(p.values,xlab="feature index",ylab="p value",main="Kolmogorov–Smirnov test p value");
b<-ggplot()+geom_point(aes(x=seq(1:length(p.values)),y=p.values))+xlab("feature index")+scale_x_continuous(breaks=seq(1:3))+ggtitle("Kolmogorov Smirnov test")+ylab("p value");
print(b,newpage = FALSE);
popViewport();


dev.off();
