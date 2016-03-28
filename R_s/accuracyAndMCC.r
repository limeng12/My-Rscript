#0 negative 1 positive
calculatingMCC<-function(standard, input){
  tp=0;  fp=0;
  fn=0;  tn=0;
  
  for(i in 1:length(standard)){
    if((standard[i]==0)&(input[i]==0) )
      tn=tn+1;
    if((standard[i]==0)&(input[i]==1) )
      fp=fp+1;
    if((standard[i]==1)&(input[i]==0) )
      fn=fn+1;
    if((standard[i]==1)&(input[i]==1) )
      tp=tp+1;
  }
  cat(str_c("tp:",tp," tn:",tn," fp:",fp," fn:",fn,"\n","auc:",(tp+tn)/(tp+tn+fp+fn),"\n" ));
  mcc<-(tp*tn-fp*fn)/(  sqrt(  (tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)  )  ) 
  
  return(mcc);
  
}

calculatingAUC<-function(standard, input){
  tp=0;  fp=0;
  fn=0;  tn=0;
  
  for(i in 1:length(standard)){
    if((standard[i]==0)&(input[i]==0) )
      tn=tn+1;
    if((standard[i]==0)&(input[i]==1) )
      fp=fp+1;
    if((standard[i]==1)&(input[i]==0) )
      fn=fn+1;
    if((standard[i]==1)&(input[i]==1) )
      tp=tp+1;
  }
  cat(str_c("tp:",tp," tn:",tn," fp:",fp," fn:",fn,"\n"));
  auc<-(tp+tn)/( tp+tn+fp+fn ) 
  return(auc);
}