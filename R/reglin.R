reglin <-
function(x, y) {
  # x and y are two vectors of length n
  n <- length(x)
  xbar <- mean(x)
  ybar <- mean(y)
  x2bar <- mean(x ^ 2)
  y2bar <- mean(y ^ 2)
  xybar <- mean(x * y)
  beta1hat <- (xybar - xbar * ybar) / (x2bar - xbar ^ 2)
  beta0hat <- ybar - beta1hat * xbar
  sigmahat <- sqrt(n * (y2bar + beta0hat ^ 2 + beta1hat ^ 2 * x2bar  
                          - 2 * beta0hat * ybar - 2 * beta1hat * xybar 
                          + 2 * beta0hat * beta1hat * xbar) / (n - 2))
  sigmahat.beta1hat <- sigmahat / sqrt(n * (x2bar - xbar ^ 2))
  test.statistic <- beta1hat / sigmahat.beta1hat
return(list(beta0hat = beta0hat, beta1hat = beta1hat, sigmahat = sigmahat,
            sigmahat.beta1hat = sigmahat.beta1hat,
            test.statistic = test.statistic))
}
