# sample
x <- rt(1000, df=4)

# Hill estimator, equation (5.23)
alpha_hat <- function(x, k) {
  x <- sort(x, decreasing = TRUE) # look at the k largest values
  return(1 / (mean(log(x[1:k])) - log(x[k])))
}


install.packages("tea")
data("danish", package = "tea")

plot(sapply(2:60, function(k) alpha_hat(danish, k)), type="l")

# To do: confidence interval based on likelihood ratio test

