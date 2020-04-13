#include "prox_fantope.h"

// min  -tr(AX) + (0.5 / alpha) * ||X - B||_F^2
// s.t. X in Fantope

//' Proximal operator on Fantope
//'
//' This function solves the optimization problem
//' \deqn{\min\quad-tr(AX) + \frac{1}{2\alpha}||X - B||_F^2}{min  -tr(AX) + (0.5 / \alpha) * ||X - B||_F^2}
//' \deqn{s.t.\quad X\in \mathcal{F}^d}{s.t. X in F^d}
//'
//' @param A       A symmetric matrix.
//' @param B       A symmetric matrix of the same dimension as \code{A}.
//' @param alpha   Proximal parameter.
//' @param d       Fantope parameter.
//' @param eps     Precision of the result.
//' @param inc     How many incremental eigenvalues to compute in each iteration.
//' @param maxiter Maximum number of iterations.
//' @param verbose Level of verbosity.
// [[Rcpp::export]]
Rcpp::NumericMatrix prox_fantope(MapMat A, MapMat B, double alpha, int d,
                                 double eps = 1e-5, int inc = 1, int maxiter = 10,
                                 int verbose = 0)
{
    const int n = A.rows();
    if(A.cols() != n)
        Rcpp::stop("A is not square");
    if(B.rows() != n || B.cols() != n)
        Rcpp::stop("dimensions of A and B do not change");

    MatrixXd mat = B + alpha * A;
    MapConstMat matm(mat.data(), n, n);

    Rcpp::NumericMatrix res(n, n);
    MapMat resm(res.begin(), n, n);
    double dsum;

    prox_fantope_hard_impl(matm, d, inc, maxiter, resm, dsum, eps, verbose);

    return res;
}
