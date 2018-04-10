gaussianize <- function (X,jitter=FALSE){ 
  #   Transform each column of data matrix X to normality using the inverse
  #   Rosenblatt transform.
  #
  # inspired by split.m in normal.m by Van Albada, S.J., Robinson P.A. (2006)
  # Transformation of arbitrary distributions to the normal distribution with application to EEG
  # test-retest reliability. J Neurosci Meth, doi:10.1016/j.jneumeth.2006.11.004
  #
  #  Written 26/06/2015 by Julien Emile-Geay (USC)
  #translated to R and added jitter option by 29/06/2015 by Nick McKay (NAU) 
  
  if(!is.matrix(X)){
    X=as.matrix(X)
  }
  p=NCOL(X)
  n=NROW(X) 
  
  if(jitter){
    #add tiny random numbers to avoid ties
    X=array(rnorm(p*n,mean=0,sd=sd(as.vector(X))/1e6),c(n,p))+X
  }
  
  
  Xn    = matrix(0,n,p);
  for (j in 1:p){
    # Sort the data in ascending order and retain permutation indices
    R=rank(X[,j])
    # The cumulative distribution function
    CDF = R/n - 1/(2*n);
    # Apply the inverse Rosenblatt transformation
    Xn[,j] = qnorm(CDF)  # Xn is now normally distributed
  }
  
  return(Xn)
}
