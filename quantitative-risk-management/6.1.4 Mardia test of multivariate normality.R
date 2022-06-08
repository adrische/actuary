# 6.1.4 Mardia test of multinormality

data(DJ, package="QRM")

stocks <- c("AXP", "EK", "BA", "C", "KO", "MSFT", "HWP", "INTC", "JPM", "DIS")

DJ <- window(DJ, start="1993-01-01", end="2000-12-29")[, stocks]

daily     <- diff(log(DJ), trim=T)
weekly    <- diff(log(daily2weekly(DJ, startOn="Mon")), trim=T)
monthly   <- diff(log(daily2monthly(DJ)), trim=T)
quarterly <- aggregate(daily, 
                       by=timeSequence(from=start(daily), to=end(daily), by="quarter"), 
                       FUN=sum)


mardia_test <- function(rets) {
  n <- nrow(rets)
  d <- ncol(rets)
  
  means <- matrix(colMeans(rets, na.rm=T), ncol=d, nrow=n, byrow = TRUE)
  
  # empirical covariance estimator (biased) (6.9)
  S   <- (1/nrow(rets))     * t(rets - means) %*% (rets - means)
  # unbiased covariance estimator, this seems identical to cov(rets)
  # S_u <- (1/(nrow(rets)-1)) * t(rets - means) %*% (rets - means)
  
  # Mahalanobis angle
  D_ij <- (rets - means) %*% solve(S) %*% t(rets - means)
  # Mahalanobis distance (6.15)
  D_i_squared <- diag(D_ij)
  
  # multivariate skewness? (6.16)
  print(b_d <- mean(D_ij^3))
  # multivariate kurtosis? (6.16)
  print(k_d <- mean(D_i_squared^2))
  
  # p-values (6.17)
  print(pchisq(n*b_d/6, df=d*(d+1)*(d+2)/6, lower.tail = FALSE))
  print(pnorm((k_d - d*(d+2))/sqrt(8*d*(d+2)/n), lower.tail = FALSE))
  
  # qqplot with theoretical chi-squared distribution 10 degrees of freedom
  # TODO: render the hard-coded "10" in xlab dynamically from d
  qqplot(qchisq(ppoints(n), df = d), D_i_squared, 
         pch=16, las=1,
         xlab=expression(paste(chi[10]^2, "-quantile")),
         ylab=expression(paste("Ordered ", D^2, " data")))
  qqline(D_i_squared, distribution = function(p) qchisq(p, df = d))
  
}

mardia_test(daily)
mardia_test(weekly)
mardia_test(monthly)
mardia_test(quarterly)

# TODO: the results do not 100% match the results in Table 6.1, 
# the p-values seem very sensitive to the values of the corresponding statistics
# and so to off-by-one errors or similar.
