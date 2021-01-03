##Simulate Frechet sample
X <- rfrechet(500, 2, 0, 1/2)
summary(X)
##Plots to visualize the data
par(mfrow = c(1,3))
meanLifePlot(X)
plot(ecdf(X), col="blue", main = "Empirical Distribution of X", xlab = "x", ylab = "Fn(x)" )
plot(X, xlab = "Frechet Extremes", ylab = "X = x", type ="b", col ="gray", 
     main = "Frechet Random Sample", pch =19)

abline(h = 100, col ="black", lwd = 2, lty =2)
par(mfrow =c(1,1))



plot(X[X<=1100], xlab = "Frechet Extremes", ylab = "X = x < 1100", type ="b", col ="gray", 
     main = "Frechet Random Sample under 1100", pch =19)
abline(h = 100, col ="black", lwd = 2, lty =2)

par(mfrow = c(1,2))
qqplot(rexp(200, 2), log(X[X>=150]), pch = 19, main = "Log Data QQ-Plot", ylab ="ln(x)", 
       xlab = "exponential quantiles")
qqnorm(X, pch = 19)
par(mfrow=c(1,1))

##Threshold search for the Frechet sample
threshs <- c(1000, 50, 200, 30, 500, 20)
K = 6
probs = rep(1/K, K)
threshSearch <- gibbsThresholdSearch(X, thresholds = threshs, probs = probs, 10000)
par(mfrow=c(2,3))
hist(threshSearch[[1]][1,-(1:5000)], main = "Threshold 1 probability", xlab = expression(beta[1]) )
hist(threshSearch[[1]][2,-(1:5000)],  main = "Threshold 2 probability", xlab = expression(beta[2]) )
hist(threshSearch[[1]][3,-(1:5000)],  main = "Threshold 3 probability", xlab = expression(beta[3]) )
hist(threshSearch[[1]][4,-(1:5000)],  main = "Threshold 4 probability", xlab = expression(beta[4]) )
hist(threshSearch[[1]][5,-(1:5000)],  main = "Threshold 5 probability", xlab = expression(beta[5]) )
hist(threshSearch[[1]][5,-(1:5000)],  main = "Threshold 6 probability", xlab = expression(beta[6]) )
par(mfrow=c(1,1))
betsMean = c( mean(threshSearch[[1]][1,-(1:500)]), mean(threshSearch[[1]][2,-(1:500)]),
              mean(threshSearch[[1]][3,-(1:500)]), mean(threshSearch[[1]][4,-(1:500)]), 
              mean(threshSearch[[1]][5,-(1:500)]),  mean(threshSearch[[1]][6,-(1:500)]) )
betsMean
sum(threshs*betsMean)

##Hypothesis tests for the Frechet sample
Xextr <- X[X >= 100]
alphaBY(Xextr, 950, 1,1)
m = 500
Xn <- max(X)
bot <- c(1:50)
top <- bot
for(i in 1:50){
  k= 10*(i-1)+1
  l = i*10
  bot[i] = k
  top[i] = l
  Xmax[i] = max(X[k:l])
}
summary(Xmax)

hypoTest(Xn, 1, 1/10)
quantile(X, probs =seq(0,1, 0.05))


Hx <- hypoChain(15000, Xmax, 10, 1/3, 2, 1 )
par(mfrow = c(1,2))
hist(Hx[-(1:7500),1], breaks = 25, main =expression(paste(
  "Posterior Distribution of ", xi)), xlab = expression(xi))
plot(ecdf(Hx[-(1:7500),1]), pch=19, main = expression(paste("Posterior Cumulative Probability ", xi)), 
     xlab = expression(xi), ylab = expression(paste("F(", xi,"| Xn)" )) )
par(mfrow = c(1,1))
plot(Hx[,1], type = "l")
plot(Hx[-(1:7500),1], type = "l")
quantile(Hx[-(1:7500),1])
##Clearly these data have a Heavy-Tail

##True values
trueA <- 1/2
trueT <- 200
trueG <- 0
##Posterior distribution simulation for the Frechet sample
MCMCx <- MHRWChain(10000,  c(1,4,0), X[X>=100], 100, c(1,1,1,1), 50) 
summary(MCMCx)

par(mfrow=c(2,3))
hist(MCMCx[-(1:5000),1], breaks = 25, main = expression(paste("Posterior Distirbution of ", alpha)),
     xlab = expression(alpha))
abline(v = mean(MCMCx[-(1:5000),1]), col = "blue", lwd = 2, lty =2)
abline(v = trueA, col = "red", lwd = 2)
legend("topright", legend = c("True Value", "Bayes Estimate"),
       col=c("red", "blue"), inset = .005, lty=1:2, cex=0.4)

hist(MCMCx[-(1:5000),2], breaks = 15,  main = expression(paste("Posterior Distirbution of ", theta)),
     xlab = expression(theta))
abline(v = mean(MCMCx[-(1:5000),2]), col = "blue", lwd = 2, lty =2)
abline(v = trueT, col = "red", lwd = 2)
legend("topleft", legend = c("True Value", "Bayes Estimate"),
       col=c("red", "blue"), inset = .05, lty=1:2, cex=0.4)

hist(MCMCx[-(1:5000),3], breaks = 20,  main = expression(paste("Posterior Distirbution of ", gamma)),
     xlab = expression(gamma))
abline(v = mean(MCMCx[-(1:5000),3]), col = "blue", lwd = 2, lty =2)
abline(v = trueG, col = "red", lwd = 2)
legend("topleft", legend = c("True Value", "Bayes Estimate"),
       col=c("red", "blue"), lty=1:2, cex=0.4)
#par (mfrow=c(1,1))
#par(mfrow=c(1,3))
plot(MCMCx[-(1:5000),1], type ="l", ylab = expression(alpha), 
     main = expression(paste("Time Series Plot of ", alpha)), xlab = "Iterarions")

plot(MCMCx[-(1:5000),2], type ="l", ylab = expression(theta), 
     main = expression(paste("Time Series Plot of ", theta)), xlab = "Iterarions")

plot(MCMCx[-(1:5000),1], type ="l", ylab = expression(gamma), 
     main = expression(paste("Time Series Plot of ", gamma)), xlab = "Iterarions")
par (mfrow=c(1,1))

##Frcehet vs Gumbel likelihood hypothesis test for extremes 
hypoTest <- function(max, a0, b0){
  alphas = rgamma(1000, a0, b0)
  priorOdds = 100
  null1 = dfrechet(max, alphas, 1, 0)
  null = mean(null1)
  alt = dgumbel(max, 0,1)
  BF = priorOdds*exp(log(alt) - log(null))
  if(BF >= 10){
    print("Reject H0, might not have a Heavy-Tail.")
  } else {
    print("Fail to Reject H0, no evidence against a Heavy-Tail")
  }
  return(BF)
} 
