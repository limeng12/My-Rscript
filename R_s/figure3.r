plotcoilasa<-function(){
  
  pdf("figure_3_random_coil_asa.pdf",width=10,height=10);
  coilDataset<-alldatasetNA[-1872,];
  rownames(coilDataset)<-paste0(sapply(str_split(rownames(coilDataset),"\\$"),"[",2),
                                sapply(str_split(rownames(coilDataset),"\\$"),"[",3),
                                sapply(str_split(rownames(coilDataset),"\\$"),"[",4))
  
  p1<-ggplot(coilDataset)+geom_point(aes(x=ss_8,y=asa_1,color=label),shape=1)+theme_classic()+
    xlab("ss_8 (min probability of the amino acid in coil)")+
    ylab("asa_1 (average ASA in translated amino acid sequence)")+
    geom_text(aes(x=ss_8,y=asa_1,label=rownames(coilDataset)),size=0.2 )+
    scale_colour_brewer(palette="Set1")#+scale_color_identity();
  
  print(p1);
  dev.off();
}
plotcoilasa();
