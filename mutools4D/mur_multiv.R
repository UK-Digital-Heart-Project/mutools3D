#' Mass Univariate Regression.

#' In mass multivariate analysis, Multivariate Analysis of Variance (MANOVA) is used to analyze patterns of heart activity across multiple regions or vertices.
#' Y represents a vector of, for instance, vertex coordinates over time, and X denotes a design matrix.
#' matrix Y (N = number of subjects, V = number of vertices in the ventricular mesh, K = number of principal components from the concatenated vertices over time)
#' matrix X  (N = number of subjects, C = number of variables + intercept term) and the numbers of the column variables to extract
#' @param X is the design matrix. Number of rows = number of subjects in the study, number of columns = number of vertices in the atlas. Numerical variable must be normalized to 0-mean and unit-standard deviation. Categorical variables must be coded using dummy coding. The first column should contain the intercept (all 1s).
#' @param Y is the imaging matrix. Number of rows = N, Number of columns = V, K= num principal components. 
#' @param extract is an array expressing which covariates in X you want to extract.
#' @keywords mur_multivariate regression
#' @export
#' @examples extract <- c(1,3) #extract the first and third covariate.
#' result <- mur_multivariate(X, Y, extract)
#' Etas <- result[,1]
#' fstatistics <- result[,2]
#' pvalues <- result[,3]
          
library(effectsize)

mur_multivariate <- function(X, Y, extract){
  
  # Y dimensions (num subjects, num vertices, num principal components)
  
  # Set the dimension names of Y
  dimnames(Y) <- list(
    paste0("obs.", seq_len(dim(Y)[1])),
    paste0("vert.", seq_len(dim(Y)[2])),
    paste0("pc.", seq_len(dim(Y)[3]))
  )
  
  stopifnot(length(dim(Y)) == 3)
  

  colnames(X) <- paste0("X.", seq_len(ncol(X)))
  
  
  nPoints <- dim(Y)[2]
  #number of points in the mesh/vertexes under study
  
  tMUR <- matrix(0, ncol=3*length(extract), nrow=nPoints)

  #beta and p obtained at each point 
  
  for(y in 1:nPoints){
    #do the regression for all the vertexes of the atlas
    #the for cycle is less time consuming than a foreach in this case
    

    
    dat <- as.data.frame(cbind(Y[, y, ], X))

    
    model_formula <- paste0("cbind(", paste0(dimnames(Y)[[3]], collapse = ","), ") ~ ",
                            paste0(colnames(X), collapse = "+"))
    results <- stats::manova(as.formula(model_formula), data = dat)
    ss <- summary(results)$stats
    # Calculate effect sizes
    effect_size <- eta_squared(results, partial=TRUE)
    
    for(iEx in 1:length(extract)){
      tMUR[y,1+(iEx-1)*3] <- effect_size$Eta2_partial[extract[iEx]]  # Eta squared
      tMUR[y,2+(iEx-1)*3] <- ss[extract[iEx], 3]  # approx F
      tMUR[y,3+(iEx-1)*3] <- ss[extract[iEx], 6]  # p-value
    }
    
  }
  
  return(tMUR)
  
}
