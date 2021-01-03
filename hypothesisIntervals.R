##test for the functions in the script, z are the tested data
iters = 1000
H = hypoChain(iters, z, -1, 1, 3, 1)
summary(H)
hist(H[-(1:iter/2),1], breaks = 15)
hist(H[,1], breaks = 15)
##Acceptance percentage
1 - mean(duplicated(H[-(1:iters/2),1]))

##Mean life plot function
meanLifePlot <- function(extr){
  Xm = max(extr)
  u = seq(0,Xm-1)
  mLP = seq(1,max(u))
  for(i in 1:length(u)){
    Xi = X[X>=u[i]]
    mLP[i] = mean(Xi-u[i])
  }
  return(plot(x = u, y = mLP, type = "l", ylab ="Mean Life Above u", main ="Mean Life Plot"))
}

##Bayes estimator for alpha based on the Pareto-Gamma conjugacy for the the tail index alpha
alphaBY <- function(extr, threshold, a0, b0){
  if(length(extr[extr>=threshold] != 0)){
    BY = (length(extr[extr>=threshold])+a0)/(sum(log(extr[extr>=threshold]/threshold))+b0)
  } else {
    BY = 0
  }
  return(BY)
}

##Verify the values fall into the area the function is defined  
supportCheck <- function(maxima, xi){
  if(abs(x1) > 1e-6) {
  b = 1+xi*maxima
  } else {
   b = rep(1, length(maxima)) 
  }
  s = length(b[b <= 0])
   if(s > 0 ){
    return(0)
  } else {
    return(1)
  }
}

##Generalized Extreme Value distribution likelihood function
likeGEV <- function(maxima, xi){
  if(supportCheck(maxima, xi) < 1) {
    loglike = -Inf
  }  else {
  if(abs(xi) > 1e-7) {
      a = sum(exp(-1/xi*log(1+xi*maxima)))
      b = (1/xi+1)*sum(log(1+xi*maxima))
      loglike = -(a+b)
    } else {
      a = sum(exp(-maxima) + maxima)
      loglike = -a
    }
  }
  return(loglike)
}

##Simulate a sample of frechet extremes
maxSample <- function(size){
  y = c(1:size)
  for(i in 1:size){
    x = rfrechet(10, scale = 1, loc =0, shape = 2)
    y[i] = max(x)
  }
  return(y)
}

##Gradient log-likelihood function
gradLogLike <- function(maxima, xi){
  if(supportCheck(maxima, xi) < 1){
    g = -Inf
  } else {
    if(abs(xi) > 1e-7){
  a = exp(-1/xi*sum(log(1+xi*maxima)))
  b = sum((log(1 + maxima*xi))*(1+a))
  c = sum((maxima/(1+xi*maxima))*(1+a))
  d = sum(maxima/(1+xi*maxima))
  g = 1/xi^2*b - (1/xi*c + d)}
    else {
      g = 0
    }
  }
  return(g)
}

##Markov Chain simulation for the posterior distrbution of xi, to carry out the hypothesis test
hypoChain <- function(iter, maxima, initial, c, L, Msd){
  chain = array(dim=c(iter + 1, 1))
  probs = array(dim=c(iter + 1, 1))
  chain[1] = initial
  for(i in 1:iter){
    momenta1 = rnorm(1, mean = 0, sd = Msd)
    props = chain[i]
    momenta = momenta1 + c*gradLogLike(maxima, props)/2
    for(j in 1:L){
      props = props + c*momenta
      momenta = momenta + c*gradLogLike(maxima, props)
    }
    momenta = momenta -c*gradLogLike(maxima, props)/2
    r = exp(likeGEV(maxima, props) + 1/2*momenta1^2 - (likeGEV(maxima, chain[i])
            + 1/2*momenta^2)
            )
    r[is.nan(r)] = 0
    probs[i] = min(1,r)
    if(runif(1) < r) {
      chain[i+1] = props
    } else {
      chain[i+1] = chain[i]
    }
    if(is.infinite(gradLogLike(maxima, chain[i+1]))){
      chain[i+1] = 0
    } else {
      chain[i+1] = chain[i+1]
    }
    
  }
  probs[iter+1] = runif(1)
  chain = cbind(chain, probs)
  return(chain)
}
