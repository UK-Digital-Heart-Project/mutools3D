#' Compute nearest neighbour matrix for TFCE. 
#'
#' Given a Tx3 matrix (T = number of triangles in the triangular mesh) containing vertices IDs of each triangle, this function computes a Ex2 matrix containing the mesh edges (E = numbr of edges).
#' To each mesh corresponds a unique NNmatrix, which is exploited by TFCE function to compute TFCE scores. In order to speed up TFCE execution time, repetitions like (A,B) and (B,A) are avoided.
#' @param triangles a Tx3 matrix (T = number of triangles in the triangular mesh) containing vertices IDs of each triangle.
#' @return A Ex2 matrix containing the mesh edges (E = numbr of edges)
#' @keywords NNmatrix 
#' @export
#' @examples NNmatrix <- computeNNmatrix(triangles)


computeNNmatrix <- function(triangles){
  
  NNmatrix <- matrix(ncol = 2)

  #for each vertex
  for (v in 1:nrow(coord)){ 
  
    tr4Point <- triangles[c(which(triangles[,1]==v),
                          which(triangles[,2]==v),
                          which(triangles[,3]==v)),]
    #find all the vertices to which v is connected
  
    NN <- unique(c(tr4Point))
    #remove duplicated
  
    NN <- NN[which(NN>v)] 
    #trick to avoid to have (a,b) and (b,a)
  
    if(length(NN)>0) NNmatrix <- rbind(NNmatrix, cbind(v,NN)) 
  }

  return(NNmatrix[-1,])

}