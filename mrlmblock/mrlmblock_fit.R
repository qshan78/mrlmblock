mrlmblock.fit = function(Adj_list, y, Om, Z, lambda = 0) {
  #' 
  #' Function to perform block penalized regression using a specified community membership
  #' with non-trivial covariance matrix
  #' 
  #' @param Adj_list a list with adjacency matrices of the same size n x n
  #' @param y a matrix of responses, note Y MUST BE CENTERED
  #' @param Om precision matrix of the noise
  #' @param intercept indicate whether intercept should be fitted
  #' @param Z membership matrix of size n by K specifying the community assignments
  #' @param lambda ridge penalty parameter
  #' 
  #' @return a list containing the fitted coefficents matrix \code{B}, 
  #' and the value of the objective function.
  #' 
  #' @references 
  #' 
  # fits B and b in the model
  # l(B,b) + lambda*||B||_F^2
  # subject to B = ZCZ^T
  
  Z = Matrix(Z)
  k = ncol(Z)
  n = ncol(Adj_list[[1]])
  m = length(Adj_list)
  lambda = m*lambda
  q = ncol(y)
  
  ZAZ = lapply(Adj_list, function(A) Matrix::crossprod(Z, Matrix::crossprod(A,Z)))
  X = as.matrix(t(sapply(ZAZ, to_vector)))  
  p = ncol(X)
  
  nK = colSums(Z)
  W = diag(to_vector(tcrossprod(nK)))
  Cvec = solve(lambda*kronecker(diag(q),W)+kronecker(Om,t(X)%*%X),as.vector(t(X)%*%(y)%*%Om))
  Cmat = matrix(Cvec,p,q)
  
  C = list()
  B = array(NA,c(n,n,q))
  for (i in 1:q) {
    C[[i]] = to_matrix(Cmat[,i],diag=T)
    B[,,i] = as.matrix(Z%*%C[[i]]%*%t(Z))
  }
  residual = y - X%*%Cmat
  return(list(B= B, C=C, residual = residual))
}

mrlmblock.fixed.Z.fixed.Om = function(Adj_list, y, Z, Om, lambda.ridge = 0.1){
  n = ncol(Adj_list[[1]])
  m = length(Adj_list)
  q = ncol(y)
  K = ncol(Z)
  
  y = scale(y,scale=F)
  b = attr(y,"scaled:center")
  
  mean.A = Reduce("+", Adj_list) / length(Adj_list)
  Adj_list = lapply(Adj_list, function(A)A-mean.A)
  
  reg.result = mrlmblock.fit(Adj_list, y, Om, Z, lambda = lambda.ridge)
  return(list(Z = Z, B = reg.result$B, C = reg.result$C, Om = Om, b = b, mean.A = mean.A, residual = reg.result$residual))
}