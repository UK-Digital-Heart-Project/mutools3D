#' Cluster-extend based thresholding method.
#'
#' Given a NxV imaging matrix Y (N = number of subjects, V = number of vertices in the ventricular mesh), a NxC model matrix X  (N = number of subjects, C = number of variables + intercept term) 
#' and the number of the column variables to extract, this function computes whether a vertex belongs to a significant cluster or not using a cluster-extend based thresholding method. 
#' The output is an array which stores as 1 the vertices that reached significance, 0 otherwise.
#' @param X is the design matrix. Number of rows = number of subjects in the study, number of columns = number of vertices in the atlas. Numerical varable must be normalized to 0-mean and unit-standard deviation. Categorical variables must be coded using dummy coding. The first column should contain the intercept (all 1s).
#' @param Y is the imaging matrix. Number of rows = N. Number of columns = V. 
#' @param extract is an array expressing which covariates in X you want to extract.
#' @param A A V-dimensional vector containing the area associated with a vertex, usually its Voronoi area.
#' @param NNmatrix  Nx2 matrix containing the mesh edges. Important: to speed up the execution please avoid repetitions like (A,B) and (B,A).
#' @param nPermutations number of permutations in the permutation test, default is 1000.
#' @param HC4m flag for triggering HC4m correction, default is FALSE.
#' @param parallel flag for triggering parallel computing, default is FALSE.
#' @param nCores flag for defining the number of cores to use, default is 1.
#' @param firsThr the cluster-forming threshold.
#' @return The output of this function contains a list of the vertices that reached significace, 0 otherwise.
#' @keywords mur cluster-extend thresholding
#' @export
#' @examples res = clusterExt(X, Y, extract, A, NNmatrix, nPermutations = 1000, HC4m = TRUE, nCores=1, thrFirst = 1)

clusterExt <- function(X, Y, extract, A, NNmatrix, nPermutations = 1000, HC4m = FALSE, parallel=FALSE, nCores=1, thrFirst = 1){
  
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
        
        #compute for each permutation the area of the larget cluster
        computed <- c(maxClusterArea(h=round(resMUR[,2],2), A=A, NNmatrix=NNmatrix, firsThr=thrFirst),
                                maxClusterArea(h=round(resMUR[,2],2), A=A, NNmatrix=NNmatrix, firsThr=-thrFirst))
        return(computed)
    }
    
  }else{
    
    for(iF in 1:nPermutations){
      Yper <-  Rz[sample(nrow(Rz)),] %*% Y
      #Y permuted for the Freedman and Lane procedure
      
      resMUR <- matrix(0, ncol=3*length(extract), nrow=ncol(Y))
      
      if(HC4m) resMUR <- murHC4m(X,Yper,extract)
      if(!HC4m) resMUR <- mur(X,Yper,extract)
      
      #compute for each permutation the area of the larget cluster
      computed <- c(maxClusterArea(h=round(resMUR[,2],2), A=A, NNmatrix=NNmatrix, firsThr=thrFirst),
                    maxClusterArea(h=round(resMUR[,2],2), A=A, NNmatrix=NNmatrix, firsThr=-thrFirst))
      
      if(iF==1) resP=computed
      else resP=rbind(resP,computed)
    }
    
  }
  
  closeAllConnections()
  
  #compute 95 percentile of the largest cluster area distribution fot negative thr
  negaP = sort(resP[,2])
  negaTHR = negaP[floor(0.95*length(negaP))]
  #compute 95 percentile of the largest cluster area distribution fot positive thr
  posP = sort(resP[,1])
  posTHR = posP[floor(0.95*length(posP))]
  
  
  #unpermuted data
  results <- matrix(0, ncol=3*length(extract), nrow=ncol(Y))
  if(HC4m) results <- murHC4m(X,Y,extract)
  if(!HC4m) results <- mur(X,Y,extract)
  
  h = round(results[,2],2)

  inde = c()
  j=0
  
  #APPLY THE METHOD ON thrFirst
  thr = thrFirst
  
  if(thr>0) origI <- which(h >= thr)
  if(thr<=0) origI <- which(h <= thr)
  #compute the list of h statistic that could be contained in a cluster
  #with forming threshold thr
  
  if(length(origI)>1){
    
    firstRowok <-  which(NNmatrix[,1] %in% origI)
    #row numbers that that have in the first column a origI value
    
    rows2Keep <- firstRowok[NNmatrix[firstRowok,2] %in% origI]
    #rows that have also a origI value int the second columns
    rm(origI)
    
    if(length(rows2Keep)>1){ 
      
      g = graph_from_edgelist(NNmatrix[rows2Keep,], directed = FALSE)
      ##compute the graph from them
      
      compo <- components(g)
      ##and extract the components
      
      memberships <- compo$membership
      ##for each vertex extract its membership
      
      nCluster <- which(compo$csize>1) 
      # cluster indexes of clusters with dimension > 1
      
      for(i in 1:length(nCluster)){
        clusterIndexes <- which(memberships == nCluster[i])
        sum(A[clusterIndexes])
        ## indexes of the vertexe of the cluster with label nCluster[i]
        if(sum(A[clusterIndexes]) > posTHR){
          if(j>0){
            inde = c(inde,clusterIndexes)
            j=j+1
          }
          if(j==0){
            inde = clusterIndexes
            j=j+1
          }
        }
      }
    }

  }

  #APPLY THE METHOD ON -thrFirst
  thr = -thrFirst
  
  if(thr>0) origI <- which(h >= thr)
  if(thr<=0) origI <- which(h <= thr)
  #compute the list of h statistic that could be contained in a cluster
  #with forming threshold thr
  
  if(length(origI)>1){
    
    firstRowok <-  which(NNmatrix[,1] %in% origI)
    #row numbers that that have in the first column a origI value
    
    rows2Keep <- firstRowok[NNmatrix[firstRowok,2] %in% origI]
    #rows that have also a origI value int the second columns
    rm(origI)
    
    if(length(rows2Keep)>1){ 
      
      g = graph_from_edgelist(NNmatrix[rows2Keep,], directed = FALSE)
      ##compute the graph from them
      
      compo <- components(g)
      ##and extract the components
      
      memberships <- compo$membership
      ##for each vertex extract its membership
      
      nCluster <- which(compo$csize>1) 
      # cluster indexes of clusters with dimension > 1
      
      for(i in 1:length(nCluster)){
        clusterIndexes <- which(memberships == nCluster[i])
        ## indexes of the vertexe of the cluster with label nCluster[i]
        if(sum(A[clusterIndexes]) > negaTHR){
          if(j>0){
            inde = c(inde,clusterIndexes)
            j=j+1
          }
          if(j==0){
            inde = clusterIndexes
            j=j+1
          }
        }
      }
    }
  }
  
  if(j>0) return(inde)
  else return(0)
}