rloglogis <- function(n, shape, scale ){
  u = runif(n)
  return(scale*(u/(1-u))^(1/shape))
}

Z <- rloglogis(200, 7/2, 2)
summary(Z)

quantile(Z, probs = seq(0,1, 0.05))
par(mfrow = c(1,3))
meanLifePlot(Z)
plot(ecdf(Z), col="blue", main = "Empirical Distribution of Z", xlab = "z", ylab = "Fn(z)" )
plot(Z, xlab = "Log-Logistic Extremes", ylab = "Z = z", type ="b", col ="gray", 
     main = "Log-Logistic Random Sample", pch =19)

abline(h = 5, col ="black", lwd = 2, lty =2)
par(mfrow =c(1,1))

par(mfrow = c(1,2))
qqplot(rexp(200, 2), log(Z[Z>=5]), pch = 19, main = "Log Data QQ-Plot", ylab ="ln(z)", 
       xlab = "exponential quantiles")
qqnorm(Z, pch = 19)
par(mfrow=c(1,1))
threshs <- c(15, 5, 6, 7, 9, 10)
K = 6
probs = rep(1/K, K)
threshSearch <- gibbsThresholdSearch(Z, thresholds = threshs, probs = probs, 10000)
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
length(Z[Z>=13])

Zextr <- Z[Z >= 6]
alphaBY(Z, 14, 1,1)
m = 500
Zn <- max(Z)
Zmax <- c(1:20)
for(i in 1:20){
  Zmax[i] = max(Z[ (10*(i-1)+1):(i*10)])
}
summary(Zmax)

hypoTest(Zn,1,1/10)

Hz <- hypoChain(15000, Zmax, 10, 1/3, 2, 1 )
par(mfrow = c(1,2))
hist(Hz[-(1:7500),1], breaks = 25, main =expression(paste(
  "Posterior Distribution of ", xi)), xlab = expression(xi))
plot(ecdf(Hz[-(1:7500),1]), pch=19, main = expression(paste("Posterior Cumulative Probability ", xi)), 
     xlab = expression(xi), ylab = expression(paste("F(", xi,"| Zn)" )) )
par(mfrow = c(1,1))
plot(Hz[,1], type = "l")
plot(Hz[-(1:7500),1], type = "l", ylab = expression(xi), main = "Time Series", 
     xlab = "Iterations")
quantile(Hz[-(1:7500),1])

##Clearly these data have a Heavy-Tail

trueA <- 7/2

MCMCz <- MHRWChain(10000,  c(5,2,0), Z[Z>=15], 15, c(1,1,1,1), 30) 
summary(MCMCz)

par(mfrow=c(2,3))
hist(MCMCz[-(1:5000),1], breaks = 25, main = expression(paste("Posterior Distirbution of ", alpha)),
     xlab = expression(alpha))
abline(v = mean(MCMCz[-(1:5000),1]), col = "blue", lwd = 2, lty =2)
abline(v = trueA, col = "red", lwd = 2)
legend("topright", legend = c("True Value", "Bayes Estimate"),
       col=c("red", "blue"), inset = .005, lty=1:2, cex=0.4)

hist(MCMCz[-(1:5000),2], breaks = 15,  main = expression(paste("Posterior Distirbution of ", theta)),
     xlab = expression(theta))

hist(MCMCz[-(1:5000),3], breaks = 20,  main = expression(paste("Posterior Distirbution of ", gamma)),
     xlab = expression(gamma))

#par (mfrow=c(1,1))
#par(mfrow=c(1,3))
plot(MCMCz[-(1:5000),1], type ="l", ylab = expression(alpha), 
     main = expression(paste("Time Series Plot of ", alpha)), xlab = "Iterarions")

plot(MCMCz[-(1:5000),2], type ="l", ylab = expression(theta), 
     main = expression(paste("Time Series Plot of ", theta)), xlab = "Iterarions")

plot(MCMCz[-(1:5000),1], type ="l", ylab = expression(gamma), 
     main = expression(paste("Time Series Plot of ", gamma)), xlab = "Iterarions")
par (mfrow=c(1,1))
