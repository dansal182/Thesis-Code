##Required package to run the function in this script
require(gtools)

##Experiment values and results
K = 5
probs = rep(1/K, K)
Y = rpareto(20, 2, 3/2)
meanLifePlot(Y)
thresh = c(8,9,10,11,13)
result = gibbsThresholdSearch(Y, thresh, probs = probs, 1000)
par(mfrow=c(2,3))
hist(result[[1]][1,-(1:500)])
hist(result[[1]][2,-(1:500)])
hist(result[[1]][3,-(1:500)])
hist(result[[1]][4,-(1:500)])
hist(result[[1]][5,-(1:500)])
par(mfrow=c(1,1))
betsMean = mean(result[[1]][,-(1:500)])
sum(thresh*betsMean)

##Gibbs Threshold search proposal done with posterior probabilities for each proposed threshold
gibbsThresholdSearch <- function(extr, thresholds, probs, iter){
  n = length(extremes)
  latent = rep(0, n)
  K = length(thresholds)
  probsMat = matrix(0, K, iter)
  alphasMat = matrix(0, K, iter)
  probsMat[,1] = probs
  alphasMat[,1] = sapply(thresholds, alphaBY, extr = extr, a0 = 1, b0 = 1)
  for(i in 2:iter){
    #probability of threshold j in {1,.., k} i.e. P(Z=k) = bj
    betas = probsMat[, i-1]
    alphas = alphasMat[, i-1]
    for(j in 1:n){
      bets = rep(0,K)
      for(t in 1:K){
        bets[t] = betas[t]*dpareto(extr[j], alphas[t], thresholds[t])
      }
      tsum = sum(bets)
      bets = bets/tsum
      latent[j] = sample.int(K, size = 1, prob = bets)
    }
    m = rep(0, K)
    a = rep(0, K)
    for(k in 1:K){
      m[k] = sum(latent==k)
      a[k] = 1 + m[k]
      extrK = subset(extr, latent==k)
      extrK = extrK[extrK>thresholds[k]]
      b[k] = 1 + sum(log(extrK/thresholds[k]))
    }
    for(l in 1:K){
    alphasMat[l,i] = rgamma(1, shape = a[l], rate = b[l])
    }
    probsMat[,i] = rdirichlet(1, m+1)
  }
  return(list(probsMat, alphasMat))
}
##Pareto density function
dpareto <- function(x, alpha, threshold){
  return(alpha/threshold*exp(-(alpha+1)*log(x/threshold)))
}
