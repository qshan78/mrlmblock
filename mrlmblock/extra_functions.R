uptri = function(A,diag = F){
  A[upper.tri(A,diag = diag)]
}

unfold = function(X){
  Xr = apply(X, 3, as.vector)
  return(Xr)
}

to_matrix = function(vec,diag=T){
  if(diag==T){
    n = (sqrt(8*length(vec)+1)-1)/2
  }else if(diag==F){
    n = (sqrt(8*length(vec)+1)+1)/2
  }
  network = matrix(0,n,n)
  network[upper.tri(network,diag)] = as.double(vec)
  network = network + t(network)
  if(diag==T){
    diag(network) = diag(network)/2
  }
  return(network)
}

to_vector = function(mat){
  mat = 2*mat
  diag(mat) = diag(mat)/2
  return(uptri(mat, diag =  T))
}