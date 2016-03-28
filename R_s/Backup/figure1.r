neutralData<-data.frame(value=rnorm(100000,2,4), label=rep("neutral",10000),stringAsFactors=FALSE);

hgmdData<-data.frame(value=rnorm(100000,20,4),label=rep("hgmd",10000),stringAsFactors=FALSE);


data<-rbind.data.frame(neutralData,hgmdData);
colnames(data)<-c("value","label")

p<-ggplot(data=data)+geom_density(aes(x=value,fill=label));

print(p);




