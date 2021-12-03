### Bayesian spatial logit model with Pólya Gamma prior ###
require(MASS)
require(BayesLogit)
require(spdep)
require(Matrix)
require(utils)
#' Spatial logit estimation using Polya-Gamma priors
#'
#' @param X An n x k matrix
#' @param Y An n x 1 vector of binary or share observations 
#' @param W An n x n non-negative matrix of row-stochastic spatial weight, with zeros on the main diagonal
#' @param MODEL Either "SAR" or "SDM" depending on the specification
#' @param beta_prior_mean An k x 1 vector of prior beta means, defaults to a vector of zeros
#' @param beta_prior_var An n x n matrix of prior beta variances, defaults to a diagonal matrix with 10^8 variance (relatively flat prior)
#' @param rho_a The beta binomial prior as on p. 142 of LeSage and Pace (2009) Introduction to Spatial Econometrics. CLC Press 
#' @param nn An n x 1 vector. Gives tge  number of responses per observation, defaults to 1
#' @param niter The total number of sampling steps
#' @param nretain The number of Gibbs sampling steps retained after convergence (must be nretain <= niter)
#' @param griddy_n The size of the grid for samplign rho, defaults to 100 
#' @param thinning Computing marginal effects is costly; we do it for only every thinning-nth draw from the posterior chain. 
#' @param rmin Minimum parameter space for the spatial autoregressive parameter rho; default -1
#' @param rmax Maximum parameter space for the spatial autoregressive parameter rho; default 1
#'
#' @return Returns a list with niter- nretainposterior draws for
#'   beta, rho, yhat, omega (the latent Polya-Gamma parameters), as well as direct, indirect, and total marginal effects 
logit_spatial = function(X,Y,W, MODEL = "SDM",
                         beta_prior_mean = matrix(0,ncol(X),1),
                         beta_prior_var = diag(ncol(X)) * 10^8,
                         rho_a = 1.01,nn = rep(1,nrow(X)),
                         niter = 1000,nretain = 500,
                         griddy_n = 100,thinning = 10,
                         rmin = -1, rmax = 1) {
  n = nrow(X)
  ndiscard = niter - nretain
  if (all(X[,1] == 1)) {
    INTERCEPT = TRUE
  } else {
    INTERCEPT = FALSE
  }
  if (MODEL == "SDM") {
    sar_k = ncol(X)
    if (INTERCEPT) {
      X = cbind(X,W %*% X[,-1])
    } else {
      X = cbind(X,W %*% X)
    }
    k = ncol(X)
  } else if (MODEL == "SAR") {
    k = ncol(X)
    sar_k = ncol(X)
  } else {
    stop("MODEL has to be either SAR or SDM")
  }
  
  beta_prob = function(rho,a) 1/beta(a,a) * ((1+rho)^(a-1) *(1-rho)^(a-1))/(2^(2*a - 1))
  
  # save the posterior draws here
  postb = matrix(0,k,nretain)
  postr = matrix(0,1,nretain)
  posty = matrix(0,n,nretain)
  postom = matrix(0,n,nretain)
  
  # pre-calculate some terms for faster draws
  beta_prior_var_inv = solve(beta_prior_var)
  kappa = Y - nn/2 
  # set-up for griddy gibbs
  griddy_n = 100
  Ais = array(0,c(n,n,griddy_n))
  AiXs = array(0,c(n,k,griddy_n))
  AiXKs = matrix(0,k,griddy_n)
  YAiXs  =matrix(0,griddy_n,k)
  logdets = matrix(NA,griddy_n,2)
  logdets[,2] = seq(rmin,rmax,length.out = griddy_n + 2)[-c(1,griddy_n+2)]
  rrhos = logdets[,2]
  cat("Pre-calculate griddy GIBBS...\n")
  pb <- txtProgressBar(min = 1,max=griddy_n,style = 3)
  for (ii in 1:griddy_n) {
    Ais[,,ii] = as.matrix(solve(.sparseDiagonal(n) - rrhos[ii] * W))
    tempA = .sparseDiagonal(n) - rrhos[ii] * W
    #logdets[ii,1] = log(det(tempA))
    Ai = solve(tempA)
    AiXs[,,ii] = as.matrix(Ai %*% X)
    
    AiXKs[,ii] = t(AiXs[,,ii]) %*% kappa
    YAiXs[ii,] = t(Y) %*% AiXs[,,ii]
    setTxtProgressBar(pb, ii)
  }
  cat("\n")
  
  
  # starting values (won't matter after sufficient draws)
  curr.rho = 0
  # start from OLS
  curr.beta = as.matrix(solve(t(X) %*% X) %*% t(X) %*% kappa)
  
  curr.A = .sparseDiagonal(n) - curr.rho * W
  curr.AiX = solve(.sparseDiagonal(n) - curr.rho * W) %*% X
  curr.mu = curr.AiX %*% curr.beta
  curr.xb = X %*% curr.beta
  
  ### Gibbs sampling
  cat("Gibbs sampling...\n")
  pb <- txtProgressBar(min = 1,max=niter,style = 3)
  for (iter in 1:niter) {
    # sample omega
    curr.om = rpg(n, nn,as.vector(curr.mu))
    yy = kappa/curr.om
    curr.Ay = curr.A %*% yy
    
    # # draw beta
    tx = curr.AiX * sqrt(curr.om)
    ty = yy * sqrt(curr.om)
    V = solve(beta_prior_var_inv + crossprod(tx) )
    b = V %*% (beta_prior_var_inv%*%beta_prior_mean +  crossprod(tx,ty)  )
    curr.beta = mvrnorm(1,b,V)
    curr.xb = X %*% curr.beta
    curr.mu = curr.AiX %*% curr.beta
    
    # Draw rho using griddy Gibbs
    mus = YAiXs %*% curr.beta
    summu = apply(AiXs,c(3),function(x) {sum(nn*log(1+exp(x %*% curr.beta)))})
    ll = mus - summu + log(beta_prob(rrhos,rho_a))
    den = ll
    y = rrhos
    adj = max(den)
    den = den - adj
    x = exp(den)
    isum = sum((y[-1] + y[-length(y)])*(x[-1]  - x[-length(x)])/2)
    z = abs(x/isum)
    den = cumsum(z)
    rnd = runif(1) * sum(z)
    ind = max(c(1,which(den <= rnd)))
    curr.rho = rrhos[ind]
    curr.A = .sparseDiagonal(n) - curr.rho * W
    curr.AiX = AiXs[,,ind]
    curr.AiXK = AiXKs[,ind]
    curr.mu = curr.AiX %*% curr.beta
    
    # we are past the burn-in, save the draws
    if (iter > ndiscard) {
      s = iter - ndiscard
      postb[,s] = as.matrix(curr.beta)
      postr[s] = curr.rho
      mu = as.matrix(exp(curr.mu))
      posty[,s] = mu / (1 + mu)
      postom[,s] = curr.om
      
      
    }
    setTxtProgressBar(pb, iter)
  }
  cat("\n")
  
  # calculate summary spatial effects
  post_calcs = seq(1,nretain,by = thinning)
  meanXs = apply(X,c(2),mean)
  diag_ind = which(c(diag(n)) == 1)
  post.direct = matrix(0,sar_k,length(post_calcs))
  post.indirect = matrix(0,sar_k,length(post_calcs))
  post.total = matrix(0,sar_k,length(post_calcs))
  
  cat("Spatial FX calculation...\n")
  pb <- txtProgressBar(min = 1,max=length(post_calcs),style = 3)
  for (iii in 1:length(post_calcs)) {
    s = post_calcs[iii]
    curr.beta = postb[,s]
    curr.rho = postr[s]
    ind = max(c(1,which(rrhos <= curr.rho)))
    
    vecAI = c(Ais[,,ind])
    meanXmu = vecAI %*% t(curr.beta[1:sar_k] * meanXs[1:sar_k])
    mu2 = vecAI %*% t(curr.beta[1:sar_k]) 
    if (MODEL == "SDM") {
      vecAIW = c(as.matrix(Ais[,,ind] %*% W))
      if (INTERCEPT) {
        meanXmu[,-1] = meanXmu[,-1] + vecAIW %*% t(curr.beta[(sar_k+1):k] * meanXs[(sar_k+1):k]) 
        mu2[,-1] = mu2[,-1] + vecAIW %*% t(curr.beta[(sar_k+1):k]) 
      } else {
        meanXmu = meanXmu + vecAIW %*% t(curr.beta[(sar_k+1):k] * meanXs[(sar_k+1):k]) 
        mu2 = mu2 + vecAIW %*% t(curr.beta[(sar_k+1):k]) 
      }
    } 
    ddd = exp(meanXmu)
    ddd = ddd / (1+ddd) *  mu2
    post.direct[,iii] = colSums(ddd[diag_ind,])/n
    post.total[,iii] = colSums(ddd) / n
    post.indirect[,iii] = post.total[,iii] - post.direct[,iii]
    setTxtProgressBar(pb, iii)
  }
  close(pb)
  
  return(list(postb = postb,
              postr = postr,
              posty = posty,
              postom = postom,
              post.direct = post.direct,
              post.indirect = post.indirect,
              post.total = post.total))
}
