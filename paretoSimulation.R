##Pareto sample simulation     
Y <- rpareto( 300, 5, 3/2)
summary(Y)

quantile(Y, probs = seq(0,1, 0.05))

##Plots to visualize the data
par(mfrow = c(1,3))
meanLifePlot(Y)
plot(ecdf(Y), col="blue", main = "Empirical Distribution of Y", xlab = "y", ylab = "Fn(y)" )
plot(Y, xlab = "Pareto Extremes", ylab = "Y = y", type ="b", col ="gray", 
     main = "Pareto Random Sample", pch =19)

abline(h = 30, col ="black", lwd = 2, lty =2)
par(mfrow =c(1,1))

par(mfrow = c(1,2))
qqplot(rexp(200, 2), log(Y[Y>=30]), pch = 19, main = "Log Data QQ-Plot", ylab ="ln(y)", 
       xlab = "exponential quantiles")
qqnorm(Y, pch = 19)
par(mfrow=c(1,1))

##Threshold search for the Pareto sample
threshs <- c(10, 20, 25, 30, 35, 40)
K = 6
probs = rep(1/K, K)
threshSearch <- gibbsThresholdSearch(Y, thresholds = threshs, probs = probs, 10000)
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

##Hypothesis tests for the Pareto sample
Yextr <- Y[Y >= 30]
alphaBY(Yextr, 32, 1,1)
m = 500
Yn <- max(Y)
Ymax <- c(1:30)
for(i in 1:30){
  Ymax[i] = max(Y[ (10*(i-1)+1):(i*10)])
}
summary(Ymax)

hypoTest(Yn,1, 10)

Hy <- hypoChain(15000, Ymax, 10, 1/3, 2, 1 )
par(mfrow = c(1,2))
hist(Hz[-(1:7500),1], breaks = 25, main =expression(paste(
  "Posterior Distribution of ", xi)), xlab = expression(xi))
plot(ecdf(Hz[-(1:7500),1]), pch=19, main = expression(paste("Posterior Cumulative Probability ", xi)), 
     xlab = expression(xi), ylab = expression(paste("F(", xi,"| Yn)" )) )
par(mfrow = c(1,1))
plot(Hy[,1], type = "l")
plot(Hy[-(1:7500),1], type = "l")
quantile(Hy[-(1:7500),1])
##Clearly these data have a Heavy-Tail

##True value of the tail index alpha
trueA <- 3/2

##Posterior distribution simulation for the Pareto sample
MCMCy <- MHRWChain(10000,  c(5,2,0), Y[Y>=32], 32, c(1,1,1,1), 50) 
summary(MCMCy)

par(mfrow=c(2,3))
hist(MCMCy[-(1:5000),1], breaks = 15, main = expression(paste("Posterior Distirbution of ", alpha)),
     xlab = expression(alpha))
abline(v = mean(MCMCy[-(1:5000),1]), col = "blue", lwd = 2, lty =2)
abline(v = trueA, col = "red", lwd = 2)
legend("topright", legend = c("True Value", "Bayes Estimate"),
       col=c("red", "blue"), inset = .005, lty=1:2, cex=0.4)

hist(MCMCy[-(1:5000),2], breaks = 15,  main = expression(paste("Posterior Distirbution of ", theta)),
     xlab = expression(theta))

hist(MCMCy[-(1:5000),3], breaks = 20,  main = expression(paste("Posterior Distirbution of ", gamma)),
     xlab = expression(gamma))

#par (mfrow=c(1,1))
#par(mfrow=c(1,3))
plot(MCMCy[-(1:5000),1], type ="l", ylab = expression(alpha), 
     main = expression(paste("Time Series Plot of ", alpha)), xlab = "Iterarions")

plot(MCMCy[-(1:5000),2], type ="l", ylab = expression(theta), 
     main = expression(paste("Time Series Plot of ", theta)), xlab = "Iterarions")

plot(MCMCy[-(1:5000),1], type ="l", ylab = expression(gamma), 
     main = expression(paste("Time Series Plot of ", gamma)), xlab = "Iterarions")
par (mfrow=c(1,1))
