reglin.pvalue <-
function(x, y, M = 10000) {
  tobs <- reglin(x, y)$test.statistic
  tM <- rep(NA, M)
  for (m in 1:M) {
    tM[m] <- reglin(sample(x), y)$test.statistic
  }
  pvalue <- mean(abs(tM) > abs(tobs))
return(pvalue)
}
