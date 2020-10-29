mrlmblock.predict = function(r, test.A.list){
  test.A.list = lapply(test.A.list, function(A)A-r$mean.A)
  ZAZ.test = lapply(test.A.list, function(A) Matrix::crossprod(r$Z, Matrix::crossprod(A,r$Z)))
  ZAZr.test = as.matrix(t(sapply(ZAZ.test, to_vector)))
  coeff = sapply(r$C, uptri,diag=T)
  y.hat = ZAZr.test%*%coeff
  y.hat = t(t(y.hat)+r$b)
  return(y.hat)
}
