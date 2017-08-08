#' Compute TFCE derived p-values using Freedman-Lane procedure.
#'
#' Given a NxV imaging matrix Y (N = number of subjects, V = number of vertices in the ventricular mesh), a NxC model matrix X  (N = number of subjects, C = number of variables + intercept term) 
#' and the numbers of the column variables to extract, this function computes for each variable specified in extract a couple of arrays: in the first array the TFCE derived p-values map on the mesh and in the second array a binary variable expressing whether that vertex has reached significance using a whole-image threshold as specified in the TFCE original paper. The output is a matrix with a number of 
#' columns equal to twice the lenght of extract and a number of rows equal to the number of vertices.
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
#' @return If verbOutput = 0 the output is a matrix composed by a number of columns twice the number of variables specified in extract. For each variable, a first array will contain the TFCE derived p-values for each mesh vertex and a second array will contain a binary variable expressing whether that vertex has reached significance using a whole-image threshold for multiple testing correction as specified in the TFCE original paper. This array will fill the output matrix column-wise.
#' If verbOutput = 1 the output is a list where the pval field contains the the pvalues computed at each vertex, TFCEmatrix field contains a V x nPermutations matrix containing the TFCE scores computed for each permutation and the tfceScores field is a V-dimensional vector containing the TFCE scores of the non-permuted data.
#' @keywords mur TFCE Freedman-Lane
#' @export
#' @examples TFCEresults = permFL(X, Y, extract, A, NNmatrix, nPermutations = 1000, HC4m = FALSE)

permFL <- function(X, Y, extract, A, NNmatrix, nPermutations = 1000, HC4m = FALSE, parallel=FALSE, nCores=1, verbOutput=0, E=0.5, H=2){
  
  set.seed(1234)
  #set seed for reproducibility
  
  Z <- X[,-extract]
  #compute Z (nuisance matrix)
  Rz <- diag(nrow(Z)) - Z %*%  solve(t(Z) %*% Z) %*% t(Z)
  
  #parallelization
  if(parallel){
    cl <- makeCluster(nCores)
    registerDoParallel(cl)
  
    resP <- foreach(iF=1:nPermutations, .packages='mutools3D', .combine=rbind)%dopar%{
        Yper <-  Rz[sample(nrow(Rz)),] %*% Y
        #Y permuted for the Freedman and Lane procedure
        
        resMUR <- matrix(0, ncol=3*length(extract), nrow=ncol(Y))
        
        if(HC4m) resMUR <- murHC4m(X,Yper,extract)
        if(!HC4m) resMUR <- mur(X,Yper,extract)
        
        computed <- matrix(0, ncol=ncol(Y), length(extract))
        
        for(iEx in 1:length(extract)){
          computed[iEx,] <- TFCE(h=round(resMUR[,2+(iEx-1)*3],2), A=A, NNmatrix=NNmatrix, E=E, H=H)
          #compute TFCE
        }
        return(computed)
    }
    
  }else{
    
    for(iF in 1:nPermutations){
      Yper <-  Rz[sample(nrow(Rz)),] %*% Y
      #Y permuted for the Freedman and Lane procedure
      
      resMUR <- matrix(0, ncol=3*length(extract), nrow=ncol(Y))
      
      if(HC4m) resMUR <- murHC4m(X,Yper,extract)
      if(!HC4m) resMUR <- mur(X,Yper,extract)
      
      computed <- matrix(0, ncol=ncol(Y), length(extract))
      
      for(iEx in 1:length(extract)){
        computed[iEx,] <- TFCE(h=round(resMUR[,2+(iEx-1)*3],2), A=A, NNmatrix=NNmatrix, E=E, H=H)
        #compute TFCE
      }
      
      if(iF==1) resP=computed
      else resP=rbind(resP,computed)
    }
    
  }
  
  closeAllConnections()
  
  significance <- matrix(0, ncol=length(extract)*2, nrow=ncol(Y))
  #TFCE derived p-values 
  
  results <- matrix(0, ncol=3*length(extract), nrow=ncol(Y))
  if(HC4m) results <- murHC4m(X,Y,extract)
  if(!HC4m) results <- mur(X,Y,extract)
  
  tfceScores <- list()
  #compute the residual matrix of Z
  
  for(iEx in 1:length(extract)){
    tfceScores[[iEx]] <-  TFCE(results[,2+(iEx-1)*3], A=A, NNmatrix=NNmatrix, E=E, H=H)
    #list of TFCE scores to analyise
    
    TFCEmatrix <- resP[seq(1,nrow(resP), by=length(extract)),]
    
    minimum = sort(apply(TFCEmatrix,1,min))
    if (length(which(minimum<0)>0)) { 
      thrMin = minimum[ceiling(0.05*nrow(TFCEmatrix))]
    } else {
      thrMin = 0
    } 
    
    maximum = sort(apply(TFCEmatrix,1,max))
    if (length(which(maximum>0)>0)) {
      thrMax = maximum[floor(0.95*nrow(TFCEmatrix))]
    }else{
      thrMax = 0
    } 
    
    for(a in 1:ncol(Y)){
      if(tfceScores[[iEx]][a]>=0){
        significance[a,1+(iEx-1)*2]  <- length(which(TFCEmatrix[,a]> tfceScores[[iEx]][a]))/nPermutations
        if(tfceScores[[iEx]][a] > thrMax) significance[a,2+(iEx-1)*2] = 1
      }
      
      if(tfceScores[[iEx]][a]<0){
        significance[a,1+(iEx-1)*2]  <- length(which(TFCEmatrix[,a]< tfceScores[[iEx]][a]))/nPermutations
        if(tfceScores[[iEx]][a] < thrMin) significance[a,2+(iEx-1)*2] = 1
      }
      
      if(significance[a,1+(iEx-1)*2]==0) significance[a,1+(iEx-1)*2] <- 1/nPermutations #minimum pvalue achievable.
    }
    
  }
  
  if(verbOutput==1) TFCEresults = list("pvalues" = significance, "TFCEmatrix" = TFCEmatrix, "tfceScores" = tfceScores)
  else TFCEresults = significance
  
  rm(significance)
  rm(TFCEmatrix)
  rm(tfceScores)
  
  return(TFCEresults)  
}