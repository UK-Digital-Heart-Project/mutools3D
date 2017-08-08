#' Compute TFCE derived p-values using a standard permutation testing procedure.
#'
#' Given a NxV imaging matrix Y (N = number of subjects, V = number of vertices in the ventricular mesh), a NxC model matrix X  (N = number of subjects, C = number of variables + intercept term) 
#' and the numbers of the column variables to extract, this function computes for each variable specified in extract the TFCE derived p-values map on the mesh. The output is a matrix with a number of 
#' columns equal to the lenght of extract and a rows equal to the number of vertices.
#' @param X is the design matrix. Number of rows = number of subjects in the study, number of columns = number of vertices in the atlas. Numerical varable must be normalized to 0-mean and unit-standard deviation. Categorical variables must be coded using dummy coding. The first column should contain the intercept (all 1s).
#' @param Y is the imaging matrix. Number of rows = N. Number of columns = V. 
#' @param extract is an array expressing which covariates in X you want to extract.
#' @param A A V-dimensional vector containing the area associated with a vertex, usually its Voronoi area.
#' @param NNmatrix  Nx2 matrix containing the mesh edges. Important: to speed up the execution please avoid repetitions like (A,B) and (B,A).
#' @param nPermutations number of permutations in the permutation test, default is 1000.
#' @param HC4m flag for triggering HC4m correction, default is FALSE.
#' @param parallel flag for triggering parallel computing, default is FALSE.
#' @param nCores flag for defining the number of cores to use, default is 1.
#' @param verbOutput flag for activating verbose output, default is 0 (off).
#' @return If verbOutput = 0 the output is a matrix containing in its rows the pvalues computed at each vertex and the number of colums referes to the variables specified in extract. If verbOutput = 1 the output is a list where the pval field contains the the pvalues computed at each vertex, TFCEmatrix field contains a V x nPermutations matrix containing the TFCE scores computed for each permutation and the tfceScores field is a V-dimensional vector containing the TFCE scores of the non-permuted data.
#' @keywords mur TFCE Freedman-Lane
#' @export
#' @examples TFCEresults = perm(X, Y, extract, A, NNmatrix, nPermutations = 1000, HC4m = FALSE)

perm <- function(X, Y, extract, A, NNmatrix, nPermutations = 1000, HC4m = FALSE, parallel=FALSE, nCores=1, verbOutput=0){
  
  set.seed(1234)
  #set seed for reproducibility
  
  #parallelization
  if(parallel){
    cl <- makeCluster(nCores)
    registerDoParallel(cl)
    
    resP <- foreach(iF=1:nPermutations, .packages='MUA3DP', .combine=rbind)%dopar%{
      Yper <-  Y[sample(1:nrow(Y)),]
      
      
      resMUR <- matrix(0, ncol=3*length(extract), nrow=ncol(Y))
      
      if(HC4m) resMUR <- murHC4m(X,Yper,extract)
      if(!HC4m) resMUR <- mur(X,Yper,extract)
      
      computed <- matrix(0, ncol=ncol(Y), length(extract))
      
      for(iEx in 1:length(extract)){
        computed[iEx,] <- TFCE(round(resMUR[,2+(iEx-1)*3],2), A, NNmatrix)
        #compute TFCE
      }
      return(computed)
    }
    
  }else{
    
    for(iF in 1:nPermutations){
      Yper <-  Y[sample(1:nrow(Y)),]
      #Y permuted for the Freedman and Lane procedure
      
      resMUR <- matrix(0, ncol=3*length(extract), nrow=ncol(Y))
      
      if(HC4m) resMUR <- murHC4m(X,Yper,extract)
      if(!HC4m) resMUR <- mur(X,Yper,extract)
      
      computed <- matrix(0, ncol=ncol(Y), length(extract))
      
      for(iEx in 1:length(extract)){
        computed[iEx,] <- TFCEsecond(round(resMUR[,2+(iEx-1)*3],2), A, NNmatrix)
        #compute TFCE
      }
      
      if(iF==1) resP=computed
      else resP=rbind(resP,computed)
    }
    
  }
  
  significance <- matrix(0, ncol=length(extract), nrow=ncol(Y))
  #TFCE derived p-values 
  
  results <- matrix(0, ncol=3*length(extract), nrow=ncol(Y))
  if(HC4m) results <- murHC4m(X,Y,extract)
  if(!HC4m) results <- mur(X,Y,extract)
  
  tfceScores <- list()
  #compute the residual matrix of Z
  
  for(iEx in 1:length(extract)){
    tfceScores[[iEx]] <-  TFCE(results[,2+(iEx-1)*3], A, NNmatrix) 
    #list of TFCE scores to analyise
    
    TFCEmatrix <- resP[seq(1,nrow(resP), by=length(extract)),]
    
    for(a in 1:ncol(Y)){
      if(tfceScores[[iEx]][a]>=0) significance[a,iEx]  <- length(which(TFCEmatrix[,a]> tfceScores[[iEx]][a]))/nPermutations
      if(tfceScores[[iEx]][a]<0) significance[a,iEx]  <- length(which(TFCEmatrix[,a]< tfceScores[[iEx]][a]))/nPermutations
      if(significance[a,iEx]==0) significance[a,iEx] = 1/nPermutations #minimum pvalue achievable.
    }
  }
  
  if(verbOutput==1) TFCEresults = list("pvalues" = significance, "TFCEmatrix" = TFCEmatrix, "tfceScores" = tfceScores)
  else TFCEresults = significance
  
  return(significance)
  
}