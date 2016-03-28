var<-data.frame(x=c(1,2,3,4),y=c(2,3,4,5))
#pdf("testFigure3.jpeg")

var2<-cbind.data.frame(var,c(rep("1",2),rep("2",2))  ,rep("red",nrow(var)) );

colnames(var2)<-c(  colnames(var),"t" ,"z"  );

p <- ggplot(var2) + geom_segment(data = var2,
                      aes(x = 0, y = 0, xend = x, yend = y,size=t),color=rep(var2$z,2),
                      arrow = arrow(length = unit(0.2, 'cm')));

p <- ggplot(var2) + geom_point(data = var2,
                                 aes( x = x, y = y,color=t,size=t)  )+scale_size_identity();




print(p);

#dev.off();
