#ifndef MARKERPEN_PROX_FANTOPE_H
#define MARKERPEN_PROX_FANTOPE_H

#include "common.h"
#include "quadprog.h"
#include "inceig_tridiag.h"
#include "walltime.h"

// min  -lambda'x + 0.5 * ||x||^2
// s.t. 0 <= xi <= 1
//      x[1] + ... + x[p] = d

inline double quadprog_sol_impl(const double* lambda, int p, int d, double* sol)
{
    MatrixXd Dmat = MatrixXd::Identity(p, p);
    VectorXd dvec(p);
    std::copy(lambda, lambda + p, dvec.data());

    int q = 2 * p + 1;
    MatrixXd Amat = MatrixXd::Zero(p, q);
    for(int i = 0; i < p; i++)
    {
        Amat(i, 0) = 1.0;
        Amat(i, i + 1) = 1.0;
        Amat(i, p + i + 1) = -1.0;
    }

    VectorXd bvec = VectorXd::Zero(q);
    bvec[0] = double(d);
    for(int i = 0; i < p; i++)
    {
        bvec[p + i + 1] = -1.0;
    }

    int meq = 1;

    VectorXd lagr(q), work(2 * p + (p * (p + 5)) / 2 + 2 * q + 1);
    Eigen::VectorXi iact(q);
    double crval;
    int nact;
    int iter[2];
    int ierr = 1;

    F77_CALL(qpgen2)
        (Dmat.data(), dvec.data(), &p, &p,
         sol, lagr.data(), &crval,
         Amat.data(), bvec.data(), &p, &q,
         &meq, iact.data(), &nact, iter,
         work.data(), &ierr);

    return crval;
}

/* // [[Rcpp::export]]
Rcpp::NumericVector quadprog_sol_r(Rcpp::NumericVector lambda, int d)
{
    int p = lambda.length();
    Rcpp::NumericVector sol(p);
    quadprog_sol_impl(lambda.begin(), p, d, sol.begin());
    return sol;
} */



// min  -tr(AX) + 0.5 * ||X||_F^2
// s.t. X in Fantope
inline int prox_fantope_impl(RefConstMat A, int d, int inc, int maxiter, RefMat res, double& dsum,
                             double eps = 1e-3, int verbose = 0)
{
    VectorXd theta(inc * maxiter + d + 1);
    IncrementalEig inceig;

    double t1 = get_wall_time();
    inceig.init(A, inc * maxiter + d + 1, d + 1);
    double t2 = get_wall_time();

    const VectorXd& evals = inceig.eigenvalues();
    double f = quadprog_sol_impl(evals.data(), inceig.num_computed(), d, theta.data());
    double theta_last = theta[inceig.num_computed() - 1];

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
        int nops = inceig.compute_next(inc);
        const VectorXd& evals = inceig.eigenvalues();
        double t2 = get_wall_time();

        double newf = quadprog_sol_impl(evals.data(), inceig.num_computed(), d, theta.data());
        theta_last = theta[inceig.num_computed() - 1];

        if(verbose > 1)
            Rcpp::Rcout << "  [prox_fantope_impl] f = " << f << ", nops = " << nops
                        << ", time_eig = " << t2 - t1 << std::endl << std::endl;

        // If f does not significantly decrease
        if(std::abs(newf - f) <= eps * std::abs(f))
            break;

        f = newf;
    }

    int pos = d;
    const int end = inceig.num_computed();
    for(; pos < end; pos++)
    {
        if(std::abs(theta[pos]) <= 1e-6)
            break;
    }

    if(verbose > 1)
    {
        const int nevals = inceig.num_computed();
        if(nevals <= 5)
        {
            Rcpp::Rcout << "  [prox_fantope_impl] evals = " << inceig.eigenvalues().head(nevals).transpose() << std::endl;
        } else {
            const int tail = std::min(5, nevals - 5);
            Rcpp::Rcout << "  [prox_fantope_impl] evals = " << inceig.eigenvalues().head(5).transpose() << " ..." << std::endl;
            Rcpp::Rcout << "                              " << inceig.eigenvalues().segment(nevals - tail, tail).transpose() << std::endl;
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
    inceig.compute_eigenvectors();
    t2 = get_wall_time();
    res.noalias() = inceig.eigenvectors().leftCols(pos) *
        theta.head(pos).asDiagonal() *
        inceig.eigenvectors().leftCols(pos).transpose();
    dsum = theta.head(d).sum();
    double t3 = get_wall_time();

    if(verbose > 1)
    {
        Rcpp::Rcout << "  [prox_fantope_impl] time_post1 = " << t2 - t1
                    << ", time_post2 = " << t3 - t2 << std::endl;
    }

    return pos;
}

/* // [[Rcpp::export]]
Rcpp::NumericMatrix prox_fantope(MapMat v, double alpha, MapMat S, int d, int inc, int maxiter)
{
    MatrixXd mat = v + alpha * S;
    MapConstMat matm(mat.data(), mat.rows(), mat.cols());

    Rcpp::NumericMatrix res(v.rows(), v.cols());
    MapMat resm(res.begin(), res.nrow(), res.ncol());
    double dsum;

    prox_fantope_impl(matm, d, inc, maxiter, resm, dsum);

    return res;
} */


#endif  // MARKERPEN_PROX_FANTOPE_H
