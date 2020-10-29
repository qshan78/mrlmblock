source("mrlmblock_fit.R")

mrlmblock.fit.Om = function(Adj_list, y, Z, Om, lambda.ridge = 0.1, lambda.omega = 0.1, MAX_STEPS){
  #' 
  #' Function to perform block penalized regression using a specified community membership
  #' with non-trivial covariance matrix
  #' 
  #' @param Adj_list a list with adjacency matrices of the same size n x n
  #' @param y a matrix of responses, note Y MUST BE CENTERED
  #' @param Om starting precision matrix of the noise
  #' @param Z membership matrix of size n by K specifying the community assignments
  #' @param lambda.ridge ridge penalty parameter
  #' @param lambda.omega glasso penalty parameter
  #' 
  #' @return a list containing the fitted coefficents matrix \code{B}, \code{Om}, 
  #' and the value of the objective function.
  #' 
  #' @references 
  #' 
  # fits B and Om in the model
  # l(B,Om) + lambda*||B||_F^2 + lambda*||Om||_1
  # subject to B = ZCZ^T
  n = ncol(Adj_list[[1]])
  m = nrow(y)
  q = ncol(y)
  K = ncol(Z)
  
  converge = F
  steps = 1
  objective = 100
  obj1 = obj3 = 0
  while (!converge&&(steps<MAX_STEPS)) {
    reg.result = mrlmblock.fit(Adj_list, y, Om, Z, lambda = lambda.ridge)
    glasso.result = glasso(s = t(reg.result$residual)%*%(reg.result$residual)/m,rho = lambda.omega,penalize.diagonal = F)
    Om = glasso.result$wi
    Om[lower.tri(Om)] = t(Om)[lower.tri(Om)]
    obj1 = sum(Om*(t(reg.result$residual)%*%reg.result$residual))/m
    obj2 = lambda.omega*sum(abs(Om))-log(det(Om))
    obj3 = lambda.ridge*sum((reg.result$B)^2)
    old.objective = objective
    objective = obj1+obj2+obj3
    if(abs((objective-old.objective)/old.objective)<1e-3){
      converge = T
    }else{
      steps = steps+1
    }
  }
  return(list(Om = Om, B= reg.result$B, C=reg.result$C, residual = reg.result$residual))
}
