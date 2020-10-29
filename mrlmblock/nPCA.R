nPCA = function(Adj_list,K=13){
  ## @param Adj_list adjacency matrix list
  ## @param K number of communities
  m = length(Adj_list)
  n = dim(Adj_list[[1]])[1]
  A = array(0, dim = c(n, n, m))
  for(i in 1:m) {
    A[, , i] = Adj_list[[i]]
  }
  node = 1:n
  
  Group = list()
  Dist = rep(0,50)
  
  for (rep in 1:50) {
    idx = kmeans(apply(A,c(1,2),mean),K,iter.max = 50, nstart=10)$cluster
    group = list()
    for (k in 1:K) {
      group[[k]] = node[idx==k]
    }
    
    for (iter in 1:200) {
      mu = array(0,c(K,K,m))
      mu2 = array(0,c(n,K,m))
      for (subi in 1:m) {
        for (i in 1:(K-1)) {
          for (j in (i+1):K) {
            mu[i,j,subi] = mean(A[idx==i,idx==j,subi])
            mu[j,i,subi] = mu[i,j,subi]
          }
        }
        for (i in 1:K) {
          mu[i,i,subi] = mean(uptri(A[idx==i,idx==i,subi],diag=T))
        }
        for (k in 1:K) {
          mu2[,k,subi] = mu[idx,k,subi]
        }
      }
      
      if(iter>2){
        old_sum_dist=sum_dist;
      }
      sum_dist = 0
      idx = rep(0,n)
      for (k in 1:K) {
        temp = group[[k]]
        for (i in 1:length(temp)) {
          vtx = temp[i]
          dist = rep(0,K)
          for (sub in 1:m) {
            temp3 = mu2[,,sub]
            temp3[node==vtx,] = 0
            temp2 = A[,node==vtx,sub]
            for (j in 1:K) {
              dist[j]= dist[j]+sum((temp2-temp3[,j])^2)
            }
          }
          idx[node==vtx] = which.min(dist)
          sum_dist = sum_dist+min(dist)
        }
      }
      if(iter==1){
        opt_group = group
        opt_sum_dist = sum_dist
      }
      if(iter>2 && sum_dist < opt_sum_dist){
        opt_group = group
        opt_sum_dist = sum_dist
      }
      for (k in 1:K) {
        group[[k]] = node[idx==k]
      }
      
      ##check here!
      flag = sum(unlist(lapply(group, function(v)(length(v)==0))))
      if(flag){
        cat("Warning: cluster with no nodes.")
        break
        }
      if(iter>2&& sum_dist==old_sum_dist) break
    }
    
    Dist[rep] = opt_sum_dist/2;
    Group[[rep]] = opt_group;
  }
  group = Group[[which.min(Dist)]]
  Z = matrix(0,n,K)
  for (k in 1:K) {
    Z[group[[k]],k] = 1
  }
  return(list(V=Z))
}