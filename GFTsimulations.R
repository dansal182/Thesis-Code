##95% qunatile computation for the Pareto distribution
varEtaPar <- function(eta, param, u){
  return((u-param[2])*(1-eta)^(-1/param[1]) + param[2] )
}
varEtaPar(0.95, param =c(2, 0), u=10)

##predictive distribution simulation for the 95% quantile
trimPost = array(dim=c(250, 3))
varP = c(1:250)
for(i in 1:250){
  top = 20*(i-1)+1 + 5000
  bot = 20*i + 5000
  trimPost[i,] = c(mean(MCMCx[(top:bot),1]), mean(MCMCx[(top:bot),2]),
                   mean(MCMCx[(top:bot),3]) )
  varP[i] = varEta(0.95, trimPost[i,])
}
summary(trimPost)
mean(varP)

Xt = rfrechet(500, 2, 0, 1/2)
Xtv = Xt[Xt>= 100]
length(Xtv)

##Predictive distribution simulation of Frechet samples to compare with the observed sample, and to simulate shocks 
postSim <- function(n){
  y = c(1:n)
  for(i in 1:n){
    x = c(1:250)
    for(j in 1:250){
      a = trimPost[j,1]
      t = trimPost[j,2]
      g = trimPost[j,3]
      x[j] = rfrechet(1, scale = t, loc = g, shape = a)
    }
    y[i] = mean(x)
  }
  return(y)
}
postPred = postSim(10000)

##Goodness of Fit test done via a uniform distribution for the simulated samples and observed values
Test = c(1:m)
for( i in 1:m){
 Test[i] = length(postPred[postPred <= Xtv[i] ])/length(postPred)
}
hist(Test, breaks = 5)
summary(Test)
quantile(postPred, probs = seq(0, 1, 0.025))

varPost = c(1:5000)

##Posterior estimate of the 95% quantile 
for(i in 1:5000){
  a = MCMCw[(5000+i), 1]
  g = MCMCw[(5000+i), 3]
  varPost[i] = varEtaPar(0.95, c(a,g), u = 10)
  
}
mean(varPost)
quantile(varPost, probs = seq(0,1, 0.05))
varEta(0.95, c(0.46, 202, -2 ))
varEta(0.95, c(0.67, 191.5, 0.4 ))

min(MCMCx[-(1:5000), 1])

