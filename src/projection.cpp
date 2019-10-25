#include <Rcpp.h>

using Rcpp::NumericVector;

// Project x = (x1, ..., xn) to {x: sum(x) = a, x >= 0}
// [[Rcpp::export]]
NumericVector proj_pos_simplex(NumericVector x, double a)
{
    const int n = x.length();

    // Sort x in decreasing order
    NumericVector u = Rcpp::clone(x);
    std::sort(u.begin(), u.end(), std::greater<double>());

    NumericVector cumu = Rcpp::cumsum(u);
    int j = n;
    double lambda = 0.0;
    for(; j > 0; j--)
    {
        lambda = (a - cumu[j - 1]) / j;
        if(u[j - 1] + lambda > 0.0)
            break;
    }

    for(int i = 0; i < n; i++)
        u[i] = std::max(x[i] + lambda, 0.0);

    return u;
}
