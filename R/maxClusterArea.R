#' Compute the largest cluster area in a statistical map.
#'
#' Given a statistical map on a 3D mesh and cluster-forming threshold firsThr this function computes the largest cluster area in the statistical map.
#' @param h A V-dimensional vector (V = number of vertices in the ventricular mesh) containing the values of a statistic at each vertex in the mesh.
#' @param A A V-dimensional vector containing the area associated with a vertex, usually its Voronoi area.
#' @param NNmatrix A Nx2 matrix containing the mesh edges. Important: to speed up the execution please avoid repetitions ofthe form (A,B) and (B,A) in the matrix.
#' @param firsThr The cluster-forming threshold.
#' @return The largest cluster area.
#' @keywords cluster-extent
#' @export
#' @examples maxArea = maxClusterArea(h, A, NNmatrix)
#' maxArea = maxClusterArea(h, A, NNmatrix, firsThr=2)

maxClusterArea <- function(h, A, NNmatrix, firsThr=1){
  #h - input statistic map - N-dimensional vector, N = number of vertices in the mesh.
  #A - area associated to each vertex of the mesh - N-dimensional vector. 
  #NNmatrix - two columns matrix storing the mesh edges, each row contain the vertex ID of the nodes at their extremes.
  #firsThr - cluster-forming threshold

  #check input data
  if(length(h)!=length(A)) stop('Lengths of h and A are different.')
  
  nPoints <- length(h)
  #the number of mesh vertexes 
  
  maxDim = 0
  #area of the maximum cluster
  
  thr <- firsThr
  #cluster threshold
    
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
        if(sum(A[clusterIndexes])>maxDim) maxDim = sum(A[clusterIndexes])
      }
 
    }
    
  }
  
  return(maxDim)
  
}