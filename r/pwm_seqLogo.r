require(ggplot2)

require(reshape2)
setwd("E:\\limeng\\RBP_network\\r");

ggLogo<-function(tvector){

 freqs<-matrix(data=tvector,byrow=TRUE,nrow=4,dimnames=list(c('A','C','G','T')))

freqdf <- as.data.frame(t(freqs))

freqdf$pos = as.numeric(as.character(rownames(freqdf)))

freqdf$height <- apply(freqdf[,c('A', 'C','G','T')], MARGIN=1,
                       FUN=function(x){2-sum(log(x^x,base=2))})

logodf <- data.frame(A=freqdf$A*freqdf$height, C=freqdf$C*freqdf$height,
                     G=freqdf$G*freqdf$height, T=freqdf$T*freqdf$height, 
                     pos=freqdf$pos)

lmf <- melt(logodf, id.var='pos')

quartz(height=3, width=8)

ggplot(data=lmf, aes(x=as.numeric(as.character(pos)), y=value))  +
    geom_bar(aes(fill=variable,order=value), position='stack', 
        stat='identity', alpha=0.5) +
    geom_text(aes(label=variable, size=value, order=value, vjust=value),
        position='stack') +
    theme_bw()

quartz.save('StackOverflow_5438474.png', type='png')
}

library("seqLogo")


library(seqLogo)
#setwd("../RBP_motif_region/HNRNPK/")

hnRNPk_exon<-read.table("../RBP_motif_region/HNRNPK/HNRNPK_exon_pwm.txt",skip=1,header=FALSE);

hnRNPk_intron_upstream<-read.table("../RBP_motif_region/HNRNPK/HNRNPK_intron_up_pwm.txt",skip=1,header=FALSE);

hnRNPk_intron_downstream<-read.table("../RBP_motif_region/HNRNPK/HNRNPK_intron_down_pwm.txt",skip=1,header=FALSE);

hnRNPk_exon_pwm<-makePWM(t(hnRNPk_exon)/1000);

hnRNPk_intron_upstream_pwm<-makePWM(t(hnRNPk_intron_upstream)/1000);

hnRNPk_intron_downstream_pwm<-makePWM(t(hnRNPk_intron_downstream)/1000);

pwm<-list();
pwm[[1]]<-t(hnRNPk_exon)/1000
pwm[[2]]<-t(hnRNPk_intron_upstream)/1000
pwm[[3]]<-t(hnRNPk_intron_downstream)/1000

region_names<-c("(A) hnRBPk binding in exon region","(B) hnRBPk binding intron upstream region","(C) hnRBPk binding intron downstrem region");


par(mfrow=c(1,3))

#pdf("Figure 6 hnRNPk_exon_pwm.pdf")
#seqLogo(hnRNPk_exon_pwm)
#dev.off();

#pdf("Figure 7 hnRNPk_intron_pwm.pdf")
#seqLogo(hnRNPk_intron_upstream_pwm)
#dev.off();


#pdf("Figure 8 hnRNPk_intron_pwm.pdf")
#seqLogo(hnRNPk_intron_downstream_pwm)
#dev.off();


mySeqLogo = seqLogo::seqLogo

bad = (sapply( body(mySeqLogo), "==", "grid.newpage()") |
        sapply( body(mySeqLogo), "==", "par(ask = FALSE)"))
body(mySeqLogo)[bad] = NULL

#norm = function(x) scale(x, center=FALSE, scale=colSums(x))

pdf("../result/Figure 6__hnRBPk_motif.pdf",width=10,height=5);

grid.newpage()
gl <- grid.layout(nrow=1, ncol=3 );
pushViewport(viewport(layout=gl)  );

for(i in 1:3){

   #pwm = norm(matrix(runif(32), nrow=4))

   pushViewport(viewport( layout.pos.row=1, layout.pos.col=i ) )
   mySeqLogo(pwm[[i]])
   grid.text(region_names[i], x=0.5, y=0.95, hjust=0.5, vjust=1);
   
   popViewport()

}

dev.off();

