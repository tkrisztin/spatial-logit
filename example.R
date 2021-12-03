rm(list=ls())

require(MASS)
require(BayesLogit)
require(spdep)
require(Matrix)
require(utils)

### generate artifical data
n=400
sar_k = 3
k = 5
X = cbind(1, matrix(rnorm(n*(sar_k-1)),n,sar_k-1))
X[,-1] = scale(X[,-1],center = TRUE, scale = FALSE)
# 7 nearest neighbour W construction from random pattern
xy <- cbind(runif(n),runif(n))
nnb =  knn2nb(knearneigh(xy, k=8))
W=Matrix(nb2mat(nnb,style='W'))
XX = cbind(X, W %*% X[,-1])

RHO = 0.5
BETA = c(0,-10,10,-5,5)
AI = as.matrix(solve(diag(n) - RHO * W))
nn = rep(1,n) # ´number of trials

MODEL = "SDM"

MU = AI %*% (XX %*% BETA + rnorm(n,0,1) ) 
pr = as.matrix(1/(1+exp(-MU)) )        # pass through an inv-logit function
Y = rbinom(n,nn,pr)      # bernoulli response variable

if (all(X[,1] == 1)) {
  INTERCEPT = TRUE
} else {INTERCEPT = FALSE}
meanXs = apply(XX,c(2),mean)
diag_ind = which(c(diag(n)) == 1)
vecAI = c(as.matrix(AI))
meanXmu = vecAI %*% t(BETA[1:sar_k] * meanXs[1:sar_k])
mu = vecAI %*% t(BETA[1:sar_k]) 
if (MODEL == "SDM") {
  vecAIW = c(as.matrix(AI %*% W))
  if (INTERCEPT) {
    meanXmu[,-1] = meanXmu[,-1] + vecAIW %*% t(BETA[(sar_k+1):k] * meanXs[(sar_k+1):k]) 
    mu[,-1] = mu[,-1] + vecAIW %*% t(BETA[(sar_k+1):k]) 
  } else {
    meanXmu = meanXmu + vecAIW %*% t(BETA[(sar_k+1):k] * meanXs[(sar_k+1):k]) 
    mu = mu + vecAIW %*% t(BETA[(sar_k+1):k]) 
  }
} 
ddd = exp(meanXmu)
ddd = ddd / (1+ddd) *  mu
DIRECT_FX = colSums(ddd[diag_ind,])/n
TOTAL_FX = colSums(ddd) / n
INDIRECT_FX = TOTAL_FX - DIRECT_FX

source("logit_spatial_sampler.R")
res = logit_spatial(X,Y,W)

### calculate posterior mean of beta and sigma
beta_mean_hat = apply(res$postb,c(1),mean)
y_mean_hat = apply(res$posty,c(1),mean)
rho_post_mean = mean(res$postr)


# spatial effects estimates
direct_post_mean = apply(res$post.direct,c(1),mean)
indirect_post_mean = apply(res$post.indirect,c(1),mean)
total_post_mean = apply(res$post.total,c(1),mean)
direct_post_sd = apply(res$post.direct,c(1),sd)
indirect_post_sd = apply(res$post.indirect,c(1),sd)
total_post_sd = apply(res$post.total,c(1),sd)
