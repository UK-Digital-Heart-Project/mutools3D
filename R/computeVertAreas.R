#' Compute 3D mesh vertices associated areas.
#'
#' Given a Tx3 matrix (T = number of triangles in the triangular mesh) containing vertices IDs of each triangle, and a Vx3 matrix  (V = number of vertices in the triangular mesh) containing vertices coordinates, this function computes the Voronoi area associated to each vertex.
#' @param triangles a Tx3 matrix (T = number of triangles in the triangular mesh) containing vertices IDs of each triangle.
#' @param coord a Vx3 matrix  (V = number of vertices in the triangular mesh) containing vertices coordinates.
#' @return A V-dimensional vector containing vertices areas.
#' @keywords voronoi 3D mesh vertex areas
#' @export
#' @examples
#' vertAreas <- computeVertAreas(triangles, coord)

computeVertAreas <- function(triangles, coord){
  
  areas <- array(0, nrow(coord))
  #array that will contain the areas, nrow = number of vertexes
  
  #iterate on all the vertexes
  for(v in 1:nrow(coord)){ 
    
    listV <- c(which(triangles[,1]==v),
               which(triangles[,2]==v),
               which(triangles[,3]==v))
    
    
    if(length(listV)>1){
      
      tr4Point <- triangles[listV,]
      
      #triangles in which the vertex is used
      #each row, a triangle - 3 columns with the vertexes ID.
      
      
      for(iTR in 1:nrow(tr4Point)){
        #iterates on all the triangles where v belongs
        
        v1 <- tr4Point[iTR,which(tr4Point[iTR,]!=v)][1]  
        #one vertex in the iTR triangle in which v belongs.
        v2 <- tr4Point[iTR,which(tr4Point[iTR,]!=v)][2]  
        #the second vertex in the iTR triangle in which v belongs
        
        pv1 <- as.matrix(coord[v1,])
        pv2 <- as.matrix(coord[v2,])
        #extract v1 and v2 coordinates
        pv <- as.matrix(coord[v,])
        #extract the coordinates of the vertex v under study
        
        #compute midpoints between v and v1 or v2
        pm1 <- c(((pv[1]+pv1[1])/2),
                 ((pv[2]+pv1[2])/2),
                 ((pv[3]+pv1[3])/2))
        
        pm2 <- c(((pv[1]+pv2[1])/2),
                 ((pv[2]+pv2[2])/2),
                 ((pv[3]+pv2[3])/2))
        
        #triangle centroid computed using the coordinates of the vertex v under study
        #and the coordinates of the other two vertexes v1 and v2 in the triangle under study
        pcentr <- c((pv1[1]+pv2[1]+pv[1])/3,(pv1[2]+pv2[2]+pv[2])/3,(pv1[3]+pv2[3]+pv[3])/3)
        
        pointsQ <- matrix(ncol = 3, nrow = 6)
        
        pointsQ[1,] <- pv  #VERTEX UNDER STUDY
        pointsQ[2,] <- pv1 #SECOND VERTEX IN THE STARTING TRIANGLE
        pointsQ[3,] <- pv2 #THIRD VERTEX IN THE STARTING TRIANGLE
        pointsQ[4,] <- pm1 #MIDDLE POINT 1
        pointsQ[5,] <- pm2 #MIDDLE POINT 2
        pointsQ[6,] <- pcentr  #TRIANGLE CENTROID
        
        #compute triangle area
        areas[v] <- areas[v] + abs(det(pointsQ[c(1,4,6),])) + abs(det(pointsQ[c(1,5,6),]))
        
      }
      
    }
    
  }
  
 return(areas)
 
}