##Process' posterior distibution simularion via proposed metropolis hastings algorithm
MC = MHRWChain(10000, c(1,2,0), X, 10, c(1,1,1,1), 25)
summary(MC)
par(mfrow=c(2,2))
hist(MC[-(1:5000),1])
hist(MC[-(1:5000),2])
hist(MC[-(1:5000),3])
par(mfrow=c(1,1))

##Likelihood function for the process
likelihoodRW <- function(param, extr, threshold, np){
  alpha = param[1]
  theta = param[2]
  u = threshold
  gam=param[3]
  m = length(extr)
  bx=sum(log((extr-gam)/theta))
  bx[is.nan(bx)] = -Inf
  loglike = m*log(alpha)-(np*((u-gam)/theta)^(-alpha)+(alpha+1)*bx)
  return(loglike)  
}
##Parameter Alpha posterior distribution
postAlphaRW <- function(param, extr, threshold, hyperParam, np){
  a0 = hyperParam[1]
  b0 = hyperParam[2]
  alphaPrior = (a0-1)*log(param[1])-param[1]*b0
  return(likelihoodRW(param,extr, threshold, np)+alphaPrior)
}
##Parameter Theta posterior distribution
postThetaRW <- function(param, extr, threshold, hyperParam, np){
  p0 = hyperParam[1]
  q0 = hyperParam[2]
  u = threshold
  thetaPrior = (p0-1)*log(param[2])-q0*param[2]
  return(likelihood(param, extr, threshold, np)+thetaPrior)
}
##Parameter Gamma posterior distribution
postGamRW <- function(param, extr, threshold, np){
  gamPrior = 0
  return(likelihoodRW(param, extr, threshold, np)+gamPrior)
}
##Alpha proposal simulation
jumpAlphaRW <- function(ax, bx){
  return(rgamma(1, shape = ax, rate = bx))
}
##Theta proposal simulation
jumpThetaRW <- function(theta, hyperParam){
  g = rnorm(1, mean = theta, sd = 1)
  if(2*(hyperParam[1]-hyperParam[2]) < g){
    g = 1 
  } else {
    g = g
  }
  return(g)
}
##Gamma proposal simulation
jumpGammaRW <- function(gam, threshold){
  t = rgumbel(1, loc = gam,  scale = 1/2)
  
  return(t)
}
##Alpha proposal density function
propsAlphaDensRW <- function(alpha, ax, bx){
  return(log(alpha)*(ax-1)-bx*alpha)
}
##Gamma proposal density
densGam <- function(gam, postr){
  return(dgumbel(gam, mu = postr[1], sig = postr[2]))
}
##Density ratio to compare the alpha proposal and current alpha in the chain
ratioAlphaRW <- function(param, extr, threshold, props, np, prior, postr){
    r = exp(postAlphaRW(param = c(props, param[2], param[3]), extr, threshold, 
                        hyperParam = prior, np) + 
              propsAlphaDensRW(param[1], postr[1], postr[2]) -
              postAlphaRW(param = param, extr, threshold, 
                          hyperParam = prior, np) -
              propsAlphaDensRW(props, postr[1], postr[2])
  )
  r[is.nan(r)] = 0
  return(min(1,r))
}

##Density ratio to compare the theta proposal and current theta in the chain
ratioThetaRW <- function(param, extr, threshold, props, np){
 
  r = exp( likelihoodRW(param = c(param[1], props, param[3]), extr, threshold
                      , np) -
             likelihoodRW(param, extr, threshold
                        , np) )
  r[is.nan(r)] = 0
  return(min(1,r))
}
##Density ratio to compare the gamma proposal and current gamma in the chain
ratioGamRW <- function(param, extr, threshold, props, np, scale){
  r = exp(
    postGamRW(param = c(param[1], param[2], props), extr, threshold, np) + 
      densGam(param[3], c(props, scale)) -
      postGamRW(param = param, extr, threshold, np) -
      densGam(props, c(param[3], scale)) - props^2 + param[3]^2
  )
  r[is.nan(r)]=0
  return(min(1,r))
}
##Random Walk to generate the parameter simulations and approximate the process' posterior distribution
MHRWChain <- function(iter, initial, extr, threshold, hyperParam, np){
  chain = array(dim=c(iter+1, 3))
  probabilities = array(dim=c(iter+1,3))
  chain[1,] = initial
  ax = hyperParam[1] + length(extr)
  bx = hyperParam[2] + sum(log(extr/threshold))
  for(i in 1:iter){
    props1 = jumpAlphaRW(ax, bx)
    r1 = ratioAlphaRW(chain[i,], extr, threshold,
          props = props1, np, c(hyperParam[1], hyperParam[2]), c(ax, bx))
    probabilities[i,1] = r1
    if(runif(1, 0, 1)< r1){
      chain[i+1,1] = props1
    } else {
      chain[i+1,1] = chain[i,1]
    }
    props2 = jumpThetaRW(chain[i,2], c(threshold, chain[i,3]))
    r2 = ratioThetaRW(param=c(chain[i+1,1], chain[i,2], chain[i,3]),
                      extr, threshold = threshold, props = props2, np)
    
    probabilities[i,2] = r2
    if(runif(1) < r2){
      chain[i+1,2] = props2
    } else {
      chain[i+1,2] = chain[i, 2]
    }
    
    props3 = jumpGammaRW(chain[i,3], threshold)
    r3 = ratioGamRW(param=c(chain[i+1,1], chain[i+1,2], chain[i,3]),
                    extr, threshold, props = props3, np, scale = 1)
    probabilities[i,3] = r3
    if(runif(1,0,1) < r3){
      chain[i+1,3] = props3
    } else {
      chain[i+1,3] = chain[i,3]
    }
    bt = sum(log((extr -chain[i+1,3])/chain[i+1,2]))
    if(is.nan(bt)){
      chain[i+1, 2] = 1
      chain[i+1, 3] = 0
    } else {
      chain[i+1,] = chain[i+1, ] 
    }
  }
  probabilities[iter+1,] = c(0,0,0)
  chain = cbind(chain, probabilities)
  return(chain)
}





