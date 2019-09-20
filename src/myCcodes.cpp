#include <RcppArmadillo.h>
using namespace Rcpp;
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
List reglinC(NumericVector x, NumericVector y) {
 int n = x.size();
 double xbar = mean(x);
 double ybar = mean(y);
 double x2bar = mean(pow(x, 2));
 double y2bar = mean(pow(y, 2));
 double xybar = mean(x * y);
 double beta1hat = (xybar - xbar * ybar) / (x2bar - pow(xbar, 2));
 double beta0hat = ybar - beta1hat * xbar;
 double sigmahat = sqrt(n * (y2bar + pow(beta0hat, 2) + pow(beta1hat, 2) * x2bar
                   - 2 * beta0hat * ybar - 2 * beta1hat * xybar
                   + 2 * beta0hat * beta1hat * xbar) / (n - 2));
 double sigmahatbeta1hat = sigmahat / sqrt(n * (x2bar - pow(xbar, 2)));
 double teststatistic = beta1hat / sigmahatbeta1hat;

return List::create(
 _["beta0hat"] = beta0hat, 
 _["beta1hat"] = beta1hat, 
 _["sigmahat"] = sigmahat,
 _["sigmahat.beta1hat"] = sigmahatbeta1hat,
 _["test.statistic"] = teststatistic);

}


// [[Rcpp::export]]
double reglinpvalueC(NumericVector x, NumericVector y, int M = 10000) {
 double tobs = reglinC(x, y)["test.statistic"];
 NumericVector tM(M);
 for (int m = 0; m < M; m++) {
   tM[m] = reglinC(RcppArmadillo::sample(x, x.size(), FALSE, NULL), y)["test.statistic"];
}
 double pvalue = mean(abs(tM) > fabs(tobs));
return(pvalue);
}
