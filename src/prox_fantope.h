#ifndef MARKERPEN_PROX_FANTOPE_H
#define MARKERPEN_PROX_FANTOPE_H

#include "common.h"
#include "inceig_tridiag.h"
#include "walltime.h"

// Given x[0] >= x[1] >= ... >= x[n-1], find c such that
// sum y[i] = sum min(1, x[i] - c) = d, 0 < d <= n
//
// sumx is cached since this function is typically called sequentially in a loop
inline double find_c_min1(const double* x, int n, double sumx, int d)
{
    // Sum of x
    // double sumx = std::accumulate(x, x + n, 0.0);

    // Each time we assume y has the form
    // 1, ..., 1, x[i]-c, ..., x[n-1]-c
    // and then compute c such that the sum is d
    // If 1 >= x[i]-c, then we have found the solution
    for(int i = 0; i < n; i++)
    {
        const double c = double(i + sumx - d) / double(n - i);
        if(x[i] - c <= 1.0)
            return c;
        sumx -= x[i];
    }
    // If we don't find c in the loop, then all elements of y are 1
    // So we let x[n-1] - c = 1
    return x[n - 1] - 1.0;
}

// Given x[0] >= x[1] >= ... >= x[n-1], find z such that
// sum z[i] = sum max(0, min(1, x[i] - c)) = d, d <= n
// This function returns ||x - lambda||^2
inline double proj_cube_simplex(const double* lambda, int p, int d, double& c, double* sol)
{
    // Let y[i] = min(1, x[i] - c)
    // Each time we assume z has the form
    // y[0], y[1], ..., y[i-1], 0, 0, ..., 0
    // and then compute c such that the sum is d
    // If y[i] <= 0, then we have found the solution of c
    double suml = std::accumulate(lambda, lambda + d, 0.0);
    c = 0.0;
    for(int i = d; i <= p; i++)
    {
        c = find_c_min1(lambda, i, suml, d);
        if(i == p || lambda[i] - c <= 0.0)
            break;
        suml += lambda[i];
    }

    double obj = 0.0;
    for(int i = 0; i < p; i++)
    {
        sol[i] = std::min(1.0, std::max(0.0, lambda[i] - c));
        obj += (sol[i] - lambda[i]) * (sol[i] - lambda[i]);
    }

    return obj;
}

// min  -tr(AX) + 0.5 * ||X||_F^2
// s.t. X in Fantope
inline int prox_fantope_hard_impl(RefConstMat A, int d, int inc, int maxiter, RefMat res, double& dsum,
                                  double eps = 1e-3, int verbose = 0)
{
    VectorXd theta(inc * maxiter + d + 1);
    IncrementalEig inceig;

    double t1 = get_wall_time();
    inceig.init(A, inc * maxiter + d + 1, d + 1, 0, 0);
    double t2 = get_wall_time();

    const VectorXd& evals = inceig.largest_eigenvalues();
    double c;
    double f = proj_cube_simplex(evals.data(), inceig.num_computed_largest(), d, c, theta.data());
    double theta_last = theta[inceig.num_computed_largest() - 1];

    if(verbose > 1)
    {
        Rcpp::Rcout << "  [prox_fantope_impl] time_init = " << t2 - t1 << std::endl;
        Rcpp::Rcout << "  [prox_fantope_impl] inc = " << inc << ", maxiter = " << maxiter << std::endl;
    }

    for(int i = 0; i < maxiter; i++)
    {
        // If theta has reached zero eigenvalues
        if(std::abs(theta_last) <= eps)
            break;

        if(verbose > 1)
            Rcpp::Rcout << "  [prox_fantope_impl] iter = " << i << std::endl;
        if(verbose > 0 && i == maxiter - 1)
            Rcpp::Rcout << "  [prox_fantope_impl] maxiter = " << maxiter << " reached!" << std::endl;

        double t1 = get_wall_time();
        int nops = inceig.compute_next_largest(inc);
        const VectorXd& evals = inceig.largest_eigenvalues();
        double t2 = get_wall_time();

        double newf = proj_cube_simplex(evals.data(), inceig.num_computed_largest(), d, c, theta.data());
        theta_last = theta[inceig.num_computed_largest() - 1];

        if(verbose > 1)
            Rcpp::Rcout << "  [prox_fantope_impl] f = " << f << ", nops = " << nops
                        << ", time_eig = " << t2 - t1 << std::endl << std::endl;

        // If f does not significantly decrease
        if(std::abs(newf - f) <= eps * std::abs(f))
            break;

        f = newf;
    }

    int pos = d;
    const int end = inceig.num_computed_largest();
    for(; pos < end; pos++)
    {
        if(std::abs(theta[pos]) <= 1e-6)
            break;
    }

    if(verbose > 1)
    {
        const int nevals = inceig.num_computed_largest();
        if(nevals <= 5)
        {
            Rcpp::Rcout << "  [prox_fantope_impl] evals = " << inceig.largest_eigenvalues().head(nevals).transpose() << std::endl;
        } else {
            const int tail = std::min(5, nevals - 5);
            Rcpp::Rcout << "  [prox_fantope_impl] evals = " << inceig.largest_eigenvalues().head(5).transpose() << " ..." << std::endl;
            Rcpp::Rcout << "                              " << inceig.largest_eigenvalues().segment(nevals - tail, tail).transpose() << std::endl;
        }

        if(pos <= 5)
        {
            Rcpp::Rcout << "  [prox_fantope_impl] theta = " << theta.head(pos).transpose() << std::endl;
        } else {
            const int tail = std::min(5, pos - 5);
            Rcpp::Rcout << "  [prox_fantope_impl] theta = " << theta.head(5).transpose() << " ..." << std::endl;
            Rcpp::Rcout << "                              " << theta.segment(pos - tail, tail).transpose() << std::endl;
        }
    }

    t1 = get_wall_time();
    inceig.compute_eigenvectors(pos, 0);
    t2 = get_wall_time();
    res.noalias() = inceig.largest_eigenvectors().leftCols(pos) *
        theta.head(pos).asDiagonal() *
        inceig.largest_eigenvectors().leftCols(pos).transpose();
    double t3 = get_wall_time();

    if(verbose > 1)
    {
        Rcpp::Rcout << "  [prox_fantope_impl] time_post1 = " << t2 - t1
                    << ", time_post2 = " << t3 - t2 << std::endl;
    }

    return pos;
}


#endif  // MARKERPEN_PROX_FANTOPE_H
