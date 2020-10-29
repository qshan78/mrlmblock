source("mrlmblock_fit.R")
source("mrlmblock_fit_Om.R")
source("extra_functions.R")

mrlmblock.ridge.greedy = function(Adj_list, y, Z0, lambda.ridge = 0.1, lambda.omega = 0.1, MAX_STEPS){
  #'
  #'Function to perform multiple response regression with block constraints by running
  #'ridge regression on block level then update community assignment using greedy search
  #'
  #' @param Adj_list a list with adjacency matrices of the same size n x n
  #' @param y a matrix of responses, note Y DO NOT NEED TO BE CENTERED
  #' @param Z0 initial membership matrix of size n by K specifying the community assignments
  #' @param lambda.ridge ridge penalty parameter
  #' @param lambda.omega glasso penalty parameter
  #' @param MAX_STEPS number of maximum iterations
  #' 
  #' @return a list containing the fitted coefficents matrix \code{B}, intercept \code{b}, 
  #' and the value of the objective function.
  n = ncol(Adj_list[[1]])
  m = length(Adj_list)
  q = ncol(y)
  K = ncol(Z0)

  y = scale(y,scale=F)
  b = attr(y,"scaled:center")
  
  mean.A = Reduce("+", Adj_list) / length(Adj_list)
  Adj_list = lapply(Adj_list, function(A)A-mean.A)

  Z = Z0
  Om = diag(q)
  
  converge = F
  steps = 1
  objective = 100

  while (!converge&&(steps<MAX_STEPS)) {
    t1 = proc.time()[3]
    reg.result = mrlmblock.fit.Om(Adj_list, y, Z, diag(q),lambda.ridge, lambda.omega, MAX_STEPS=10)
    Om = reg.result$Om
    obj1 = sum(Om*(t(reg.result$residual)%*%reg.result$residual))/m
    obj2 = lambda.omega*sum(abs(Om))-log(det(Om))
    obj3 = lambda.ridge*sum((reg.result$B)^2)
    obj = obj1+obj2+obj3
    cat(paste0("Step ",steps,"\n",
               "After B/Om step: ",
               "Obj: ",obj,
               " .\n"))
    t2 = proc.time()[3]
    cat(paste0("Time taken in step 1:",t2-t1,"\n"))

    Z.result = Z.update(reg.result,Adj_list,Z,q,Om,lambda.ridge,y)
    Z.new = Z.result$Z.new
    Z.change = sum((tcrossprod(Z)-tcrossprod(Z.new))^2)
    t3 = proc.time()[3]
    cat(paste0("Time taken in step 2:",t3-t2,"\n"))
    
    obj2 = lambda.omega*sum(abs(Om))-log(det(Om))
    B = Z.result$B
    obj3 = lambda.ridge*sum((B)^2)
    res = Z.result$residual
    obj1 = sum(Om*(t(res)%*%res))/m
    obj = obj1+obj2+obj3
    cat(paste0("After Z step: ",
               " Obj: ",obj,
               ". Change in Z: ",
               Z.change,
               " .\n"))

    if(Z.change==0){
      converge = T
    }else{
      Z = Z.new
      steps = steps+1
    }
  }
  return(list(Z = Z.new, B = reg.result$B, C = reg.result$C, Om = Om, b = b, mean.A = mean.A, steps = steps, residual = Z.result$residual))
}

Z.update = function(reg.result,Adj_list,Z,q,Om,lambda.ridge,y){
  m = length(Adj_list)
  n = dim(Adj_list[[1]])[1]
  K = ncol(Z)
  nodes = sample(1:n,n)
  residual = reg.result$residual
  B = reg.result$B
  for (vtx in nodes) {
    new.objective = rep(0,K)
    current.label = which(Z[vtx,]==1)
    X = t(sapply(Adj_list, function(A) A[vtx,]))
    resid.change = list()
    B.temp = list()
    ZC = lapply(reg.result$C, function(Ci)Z%*%Ci)
    for (k in 1:K) {
      resid.change[[k]] = matrix(0,m,q)
      B.temp[[k]] = B
      for (j in 1:q) {
        B.change = Matrix(0,n,n,sparse = T)
        B.change[vtx,] = ZC[[j]][,k]-ZC[[j]][,current.label]
        B.change = B.change + t(B.change)
        B.change[vtx,vtx] = B.change[vtx,vtx]+reg.result$C[[j]][k,k]+reg.result$C[[j]][current.label,current.label]-2*reg.result$C[[j]][current.label,k]
        temp = -2*as.matrix((X%*%B.change[,vtx]))
        resid.change[[k]][,j] = temp[,1]
        B.temp[[k]][,,j] = B.temp[[k]][,,j]+as.matrix(B.change)
      }
      obj3 = lambda.ridge*sum((B.temp[[k]])^2)
      res = residual+resid.change[[k]]
      obj1 = sum(Om*(t(res)%*%res))/m
      new.objective[k] = obj1+obj3
    }
    new.label = which.min(new.objective)
    if(new.label!=current.label){
      Z[vtx,current.label]=0
      Z[vtx,new.label]=1
      residual = residual+resid.change[[new.label]]
      B = B.temp[[new.label]]
    }
    
  }
  return(list(Z.new=Z,residual=residual,B=B))
}
