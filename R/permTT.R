#' Compute TFCE derived p-values for t-tests using a standard permutation testing procedure.
#'
#' Given X1, X2 matrices (ncol=# of vertices, nrow=# number of samples) this function computes at each vertex of the atlas the TFCE-derived pvalues of the t-test at each vertex between one column of X1 and one column of X2.
#' @param X1 Input matrix 1 (ncol=# of vertices, nrow=# number of samples).
#' @param X2 Input matrix 2 (ncol=# of vertices, nrow=# number of samples).
#' @param A A V-dimensional vector containing the area associated with a vertex, usually its Voronoi area.
#' @param NNmatrix  Nx2 matrix containing the mesh edges. Important: to speed up the execution please avoid repetitions like (A,B) and (B,A).
#' @param nPermutations number of permutations in the permutation test, default is 1000.
#' @param parallel flag for triggering parallel computing, default is FALSE.
#' @param nCores flag for defining the number of cores to use, default is 1.
#' @return An array with a number of entries equal to the number of vertices under exam.
#' @keywords t-test TFCE 
#' @export
#' @examples TFCEresults = permTT(X1, X2, A, NNmatrix, nPermutations = 1000, parallel=FALSE, nCores=1)

permTT <- function(X1, X2, A, NNmatrix, nPermutations = 1000, parallel=FALSE, nCores=1){
  
  set.seed(1234)
  #set seed for reproducibility
  
  XtExt = rbind(X1,X2)
  
  #parallelization
  if(parallel){
    cl <- makeCluster(nCores)
    registerDoParallel(cl)
    
    resP <- foreach(iF=1:nPermutations, .packages='mutools3D', .combine=rbind)%dopar%{
      
      extract = sample(1:(nrow(XtExt)))
      Xper1 <-  XtExt[extract[1:nrow(X1)],]
      Xper2 <-  XtExt[extract[(nrow(X1)+1):length(extract)],]
      #inter-class permutation
      
      h_t <- sapply(1:ncol(Xper1), function(i) t.test(x=Xper1[,i],y=Xper2[,i],alternative="two.sided",mu=0,paired=FALSE,var.equal=FALSE,conf.level=0.95)$statistic )
      #perform paired t-test
      
      computed <- TFCE(round(h_t,2), A, NNmatrix)

      return(computed)
    }
    
  }else{
    
    for(iF in 1:nPermutations){
      extract = sample(1:(nrow(XtExt)))
      Xper1 <-  XtExt[extract[1:nrow(X1)],]
      Xper2 <-  XtExt[extract[(nrow(X1)+1):length(extract)],]
      #inter-class permutation
      
      h_t <- sapply(1:ncol(Xper1), function(i) t.test(x=Xper1[,i],y=Xper2[,i],alternative="two.sided",mu=0,paired=FALSE,var.equal=FALSE,conf.level=0.95)$statistic )
      #perform pairred t-test
      
      computed <- TFCE(round(h_t,2), A, NNmatrix)
      #tfce
      
      if(iF==1) resP=computed
      else resP=rbind(resP,computed)
    }
    
  }
  
  significance <- array(0, dim=ncol(X1))
  #TFCE derived p-values 
  
  results <- sapply(1:ncol(X1), function(i) t.test(x=X1[,i],y=X2[,i],alternative="two.sided",mu=0,paired=FALSE,var.equal=FALSE,conf.level=0.95)$statistic)
  
  tfceScores <-  TFCE(round(results,2), A, NNmatrix) 
  #list of TFCE scores to analyise
    
  for(a in 1:length(tfceScores)){
      if(tfceScores[a]>=0) significance[a]  <- length(which(resP[,a] > tfceScores[a]))/nPermutations
      if(tfceScores[a]<0) significance[a]  <- length(which(resP[,a] < tfceScores[a]))/nPermutations
      if(significance[a]==0) significance[a] = 1/nPermutations #minimum pvalue achievable.
    }
  
  
  return(significance)
  
}