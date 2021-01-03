##Gamma-Pareto tail mixture sample simulation

rmix <- function(n, a0, b0, u0, alpha0){
  p=pgamma(u0, rate = b0, shape = a0)
  w <- c(1:n)
  for(i in 1:n){
    if(runif(1) < p){
      w[i] = rgamma(1, shape = a0, rate = b0)
    } else {
      w[i] = rpareto(1, u= u0, alpha = alpha0)
    }
    
  }
  return(w)
}
W <- rmix(500, 3, 1, 3, 2)
summary(W)

##Plots to visualize the data
par(mfrow = c(1,3))
meanLifePlot(W)
plot(ecdf(W), col="blue", main = "Empirical Distribution of W", xlab = "w", ylab = "Fn(w)" )
plot(W, xlab = "Mixture Extremes", ylab = "W = w", type ="b", col ="gray", 
     main = "Mixture Random Sample", pch =19)

abline(h = 11, col ="black", lwd = 2, lty =2)
par(mfrow =c(1,1))

par(mfrow = c(1,2))
qqplot(rexp(200, 2), log(W[W>=10]), pch = 19, main = "Log Data QQ-Plot", ylab ="ln(w)", 
       xlab = "exponential quantiles")
qqnorm(W, pch = 19)
par(mfrow=c(1,1))

##Threshold search for the mixture simulation
threshs <- c(18, 4, 3, 5, 2, 11)
K = 6
probs = rep(1/K, K)
threshSearch <- gibbsThresholdSearch(W, thresholds = threshs, probs = probs, 10000)
par(mfrow=c(2,3))
hist(threshSearch[[1]][1,-(1:5000)], main = "Threshold 1 probability", xlab = expression(beta[1]), 
     xlim = c(0,1))
hist(threshSearch[[1]][2,-(1:5000)],  main = "Threshold 2 probability", xlab = expression(beta[2]),
     xlim = c(0,1))
hist(threshSearch[[1]][3,-(1:5000)],  main = "Threshold 3 probability", xlab = expression(beta[3]),
     xlim = c(0,1))
hist(threshSearch[[1]][4,-(1:5000)],  main = "Threshold 4 probability", xlab = expression(beta[4]),
     xlim = c(0,1))
hist(threshSearch[[1]][5,-(1:5000)],  main = "Threshold 5 probability", xlab = expression(beta[5]),
     xlim = c(0,1))
hist(threshSearch[[1]][5,-(1:5000)],  main = "Threshold 6 probability", xlab = expression(beta[6]),
     xlim = c(0,1) )
par(mfrow=c(1,1))
betsMean = c( mean(threshSearch[[1]][1,-(1:5000)]), mean(threshSearch[[1]][2,-(1:5000)]),
              mean(threshSearch[[1]][3,-(1:5000)]), mean(threshSearch[[1]][4,-(1:5000)]), 
              mean(threshSearch[[1]][5,-(1:5000)]),  mean(threshSearch[[1]][6,-(1:5000)]) )
betsMean
sum(threshs*betsMean)
length(W[W>=16])

##Hypothesis tests for the mixture simulation
Wextr <- W[W >=15]
alphaBY(W, 15.5, 1,1)
m = 500
Wn <- max(W)
Wmax <- c(1:50)
for(i in 1:50){
  Wmax[i] = max(W[ (10*(i-1)+1):(i*10)])
}
summary(Wmax)

hypoTest(Wn,1,1/10)

Hw <- hypoChain(15000, Wmax, 5, 1, 3, 2 )
par(mfrow = c(1,2))
hist(Hw[-(1:7500),1], breaks = 25, main =expression(paste(
  "Posterior Distribution of ", xi)), xlab = expression(xi))
plot(ecdf(Hw[-(1:7500),1]), pch=19, main = expression(paste("Posterior Cumulative Probability ", xi)), 
     xlab = expression(xi), ylab = expression(paste("F(", xi,"| Wn)" )) )
par(mfrow = c(1,1))
plot(Hw[,1], type = "l")
plot(Hw[-(1:7000),1], type = "l")
quantile(Hw[-(1:7000),1], probs = seq(0,1, 0.05))
##Clearly these data have a Heavy-Tail

##Mixture simulation tail index value
trueA <- 2

##Posterior distribution approximation 
MCMCw <- MHRWChain(10000,  c(5,2,0), W[W>=11], 11, c(1,1,1,1), 50) 
summary(MCMCw)

par(mfrow=c(2,3))
hist(MCMCw[-(1:5000),1], breaks = 25, main = expression(paste("Posterior Distirbution of ", alpha)),
     xlab = expression(alpha))
abline(v = mean(MCMCw[-(1:5000),1]), col = "blue", lwd = 2, lty =2)
abline(v = trueA, col = "red", lwd = 2)
legend("topright", legend = c("True Value", "Bayes Estimate"),
       col=c("red", "blue"), inset = .005, lty=1:2, cex=0.4)

hist(MCMCw[-(1:5000),2], breaks = 15,  main = expression(paste("Posterior Distirbution of ", theta)),
     xlab = expression(theta))

hist(MCMCw[-(1:5000),3], breaks = 20,  main = expression(paste("Posterior Distirbution of ", gamma)),
     xlab = expression(gamma))

#par (mfrow=c(1,1))
#par(mfrow=c(1,3))
plot(MCMCw[-(1:5000),1], type ="l", ylab = expression(alpha), 
     main = expression(paste("Time Series Plot of ", alpha)), xlab = "Iterarions")

plot(MCMCw[-(1:5000),2], type ="l", ylab = expression(theta), 
     main = expression(paste("Time Series Plot of ", theta)), xlab = "Iterarions")

plot(MCMCw[-(1:5000),1], type ="l", ylab = expression(gamma), 
     main = expression(paste("Time Series Plot of ", gamma)), xlab = "Iterarions")
par (mfrow=c(1,1))
