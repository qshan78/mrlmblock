library(MASS)
library(Matrix)
library(glasso)
library(RSpectra)

source("extra_functions.R")
source("mrlmblock_predict.R")
source("mrlmblock_ridge_greedy.R")
source("nPCA.R")

## generate data
y = read.table("pheno_resid.txt", quote="\"", comment.char="")
y = unname(as.matrix(y))
A = read.csv("Amat.txt",header = FALSE, nrows=5937, colClasses = rep("numeric",87153))
A = unname(as.matrix(A))

Adj_list = list()
for (i in 1:5937) {
  Adj_list[[i]] = to_matrix(A[i,],diag = F)
}

m = 5937
n.train = 500
samp = sample(m,n.train)
tasks = 1:11


### algorithm starts here
K=13
MAX_STEPS = 2
lambda.ridge = 1e-5
lambda.omega=0.1


q = length(tasks)
train.A = A[samp,]
train.y = y[samp,tasks]
train.A.list = Adj_list[samp]

test.y = y[-samp,tasks]
test.A.list = Adj_list[-samp]


# Initialize the community assignment.
# Here I run Yura's algorithm on the marginal correlation between each
# task and edge, this can be substitute with any known initialization
# (e.g. Gordon) or running Yura's algorithm on things like neurosignatures,
# regression coefficients, etc.

sigma_AY_list = list()
for (j in 1:q) {
  sigma_AY = apply(train.A, 2, function(a) cov(a,train.y[,j]))
  sigma_AY_list[[j]] = to_matrix(sigma_AY,diag = F)
}

Z0 = nPCA(sigma_AY_list,K)$V

# main algorithm
r = mrlmblock.ridge.greedy(train.A.list, train.y, Z0, lambda.ridge, lambda.omega, MAX_STEPS)

# calculate correlation between y and yhat
hat.y = mrlmblock.predict(r,test.A.list)
hat.y.train = mrlmblock.predict(r,train.A.list)

rho.test = rep(0,q)
rho.train = rep(0,q)
for (i in 1:q) {
  rho.test[i] = cor(hat.y[,i],test.y[,i])
  rho.train[i] = cor(hat.y.train[,i],train.y[,i])
}

result = list(rho.test = rho.test,
              rho.train = rho.train,
              Z.mrlm = r$Z,
              Om.mrlm = r$Om,
              C = r$C)

saveRDS(result,file = "results.rds")