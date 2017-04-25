## ------------------------------------------------------------------------
library(mutools3D)
data(Xtest)
data(Ytest)

## ------------------------------------------------------------------------
#extract results for synthetic variable
extract = 6

#it is also possible to study more than one variable at the same time
#extract = c(5,6)

result = mur(X,Y, extract)
#or
result = murHC4m(X,Y, extract)

## ------------------------------------------------------------------------
#load data for TFCE
data(NNmatrix)
data(areas)

## ------------------------------------------------------------------------
#run TFCE on the t-statistic map previously obtained
TFCEresults = TFCE(result[,2], A, NNmatrix)

