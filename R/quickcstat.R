# Approximation of the C-statistic (Large n)

quickcstat <- function(y, pred, seed=1){
  #set.seed(seed)
  casepred=pred[y == 1]
  conpred=pred[y == 0]
  
  if (length(conpred)>length(casepred)){
       conpred=conpred[sample(length(conpred),length(casepred),replace=FALSE)]
       auc.true=sum(casepred>conpred)/length(casepred)} else 
         {
      casepred=casepred[sample(length(casepred),length(conpred),replace=FALSE)]
      auc.true=sum(casepred>conpred)/length(conpred)
         }
  
  return(auc.true)
}