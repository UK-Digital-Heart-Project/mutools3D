#' Mass Univariate Regression using HC4m estimator.
#'
#' Fit a linear regression model at each vertex of a 3D atlas. HC4m estimator is used to correct for heteroskedastic data. Inputs of the function are a NxC matrix X modeling an effect under study (N = number of subjects, C = number of variables + intercept term), 
#' a NxV imaging matrix Y containing the values of a 3D phenotype at each atlas vertex (V = number of vertices in the 3D mesh), and an array (extract) containg the positions of
#' variables in X of which extract the informations of interest. The output is a Vx(3xlength(extract)) matrix containing the regression coefficient, its related t-statistic and the p-value at each vertex of the computational model for each variable specifiec in extract.
#' @param X is the design matrix. Number of rows = number of subjects in the study, number of columns = number of vertices in the atlas. Numerical varable must be normalized to 0-mean and unit-standard deviation. Categorical variables must be coded using dummy coding. The first column should contain the intercept (all 1s).
#' @param Y is the imaging matrix. Number of rows = N. Number of columns = V. 
#' @param extract is an array expressing which covariates in X you want to extract.
#' @keywords mur regression
#' @export
#' @examples extract <- c(1,3) #extract the first and third covariate.
#' result <- murHC4m(X, Y, extract)
#' betas <- result[,1]
#' tstatistics <- result[,2]
#' pvalues <- result[,3]

murHC4m <- function(X, Y, extract){
  
  nPoints <- ncol(Y)
  #number of points in the mesh/vertexes under study
  
  n <- nrow(X)
  p <- ncol(X)
  np <- n/p
  
  tMUR <- matrix(0, ncol=3*length(extract), nrow=nPoints)
  #t-statistics vector
  
  tX <- t(X)
  #transpose X
  xTxInv <- solve(crossprod(X))
  #(tX X)^-1
  
  xii <- apply(X, 1, function(x,xTxInv) t(x) %*% xTxInv %*% x, xTxInv)
  delta = sapply(xii, function(x) min(1, np*x ) +  min(1.5,np*x))
  den=(1-xii)^delta
  
  for(y in 1:nPoints){
    #do the regression for all the vertexes of the atlas
    #the for cycle is less time consuming than a foreach in this cas
    
    beta <- xTxInv %*% tX %*% Y[,y]
    #regression coefficients
    resid <-  Y[,y] - (X %*% beta)
    #residuals
    
    for(iEx in 1:length(extract)){
      hc <- ((xTxInv[extract[iEx],] %*% tX) * t((resid*resid)/den)) %*% X %*% xTxInv[,extract[iEx]] 
      se <- sqrt(hc)
      
      t <- beta[extract[iEx]]/se
      
      tMUR[y,1+(iEx-1)*3] <- beta[extract[iEx]] #beta
      tMUR[y,2+(iEx-1)*3] <- t #t
      tMUR[y,3+(iEx-1)*3] <- 2 * pt(-abs(t), df=nrow(Y)-ncol(X)) # p-value
    }
    
  }
  
  return(tMUR)
    
}