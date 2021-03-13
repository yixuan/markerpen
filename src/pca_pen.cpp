#include "common.h"
#include "prox_fantope.h"
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>

using Rcpp::IntegerVector;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using Rcpp::List;

// For two orthogonal matrices U and V, U'U = V'V = I_d,
// ||UU' - VV'||^2 = 2 * d - 2 * ||U'V||^2
inline double projection_diff(const Matrix& u, const Matrix& v)
{
    const int d = u.cols();
    Matrix uv(d, d);
    uv.noalias() = u.transpose() * v;
    return 2.0 * d - 2.0 * uv.squaredNorm();
}

// Proximal-proximal-gradient optimizer
class PPGOptimizer
{
private:
    const int   m_p;      // Dimension of the matrix, p
    const int   m_pp;     // p^2
    RefConstMat m_linear; // Linear term in the objective function

    Matrix m_z1;          // Auxiliary variable
    Matrix m_z2;          // Auxiliary variable
    Matrix m_work;        // Work space
    Matrix m_evecs;       // Eigenvectors
    Matrix m_ework;       // Work space for computing eigenvectors

    int m_fan_inc;        // Parameter for computing the Fantope proximal operator
    int m_fan_maxinc;     // Parameter for computing the Fantope proximal operator
    int m_fan_maxiter;    // Parameter for computing the Fantope proximal operator

public:
    PPGOptimizer(const RefConstMat& linear) :
        m_p(linear.rows()), m_pp(m_p * m_p), m_linear(linear),
        m_z1(m_p, m_p), m_z2(m_p, m_p), m_work(m_p, m_p),
        m_evecs(m_p, 1)
    {}

    inline void init(const RefConstMat& x0, int fan_maxinc, int fan_maxiter, bool compute_diff_evec = true)
    {
        m_z1.noalias() = x0;
        m_z2.noalias() = m_z1;

        m_fan_inc = 2;
        m_fan_maxinc = fan_maxinc;
        m_fan_maxiter = fan_maxiter;

        if(compute_diff_evec)
            m_ework.resize(m_p, 1);
    }

    // z1 <- (z1 - z2) / 2 + prox_fantope(z2 + alpha * S)
    // Returns ||newz1 - z1||_F
    inline double update_z1(double lr, double eps, int verbose)
    {
        // Compute prox_fantope(z2 + alpha * linear)
        m_work.noalias() = m_z2 + lr * m_linear;
        m_fan_inc = prox_fantope_hard_impl(
            m_work, 1, m_fan_inc, m_fan_maxiter, m_work, eps, verbose
        );

        // Update z1
        double* z1p = m_z1.data();
        const double* z2p = m_z2.data();
        const double* pfp = m_work.data();
        double diff = 0.0;

#if defined(_OPENMP)
        #pragma omp simd aligned(z1p, z2p, pfp: 32)
#endif
        for(int i = 0; i < m_pp; i++)
        {
            const double newz1 = 0.5 * (z1p[i] - z2p[i]) + pfp[i];
            diff += std::pow(newz1 - z1p[i], 2);
            z1p[i] = newz1;
        }

        // Adjust algorithm parameter
        if(verbose > 1)
            Rcpp::Rcout << "fan_dim = " << m_fan_inc << std::endl;
        m_fan_inc = std::max(5, int(1.5 * m_fan_inc));
        m_fan_inc = std::min(m_fan_inc, m_fan_maxinc);
        m_fan_inc = std::min(m_fan_inc, int(m_p / 10));

        return std::sqrt(diff);
    }

    // z2 <- (z2 - z1) / 2 + pmax(z1, 0)
    //     = { (z1 + z2) / 2, if z1 >= 0
    //       { (z2 - z1) / 2, otherwise
    // Returns ||newz2 - z2||_F
    inline double update_z2(double lr, int verbose)
    {
        const double* z1p = m_z1.data();
        double* z2p = m_z2.data();
        double diff = 0.0;

#if defined(_OPENMP)
        #pragma omp simd aligned(z1p, z2p: 32)
#endif
        for(int i = 0; i < m_pp; i++)
        {
            const double z1 = z1p[i];
            const double z2 = z2p[i];
            const double newz2 = (z1 > 0.0) ? (0.5 * (z1 + z2)) : (0.5 * (z2 - z1));
            diff += std::pow(newz2 - z2, 2);
            z2p[i] = newz2;
        }
        return std::sqrt(diff);
    }

    // Returns ||UU' - VV'||_F, U = new eigenvectors, V = old eigenvectors
    // m_ework needs to be initialized by init(compute_diff_evec = true)
    inline double update_evecs(bool first_iter = false)
    {
        Matrix& x = m_work;
        x.noalias() = 0.5 * (m_z1 + m_z2);
        Spectra::DenseSymMatProd<double> op(x);
        Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double> >
            eigs(&op, 1, 10);
        eigs.init();
        eigs.compute();

        m_ework.noalias() = eigs.eigenvectors();
        double diff = std::numeric_limits<double>::infinity();
        // Only compute the difference when ther exist old eigenvectors
        if(!first_iter)
            diff = projection_diff(m_ework, m_evecs);
        m_ework.swap(m_evecs);
        return std::sqrt(diff);
    }

    // On convergence, x = pmax(z1, 0)
    const Matrix& get_saprse_x(double lr, bool update_ev = true)
    {
        Matrix& x = m_work;
        x.noalias() = m_z1.cwiseMax(0.0);

        if(update_ev)
        {
            Spectra::DenseSymMatProd<double> op(x);
            Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double> >
                eigs(&op, 1, 10);
            eigs.init();
            eigs.compute();
            m_evecs.noalias() = eigs.eigenvectors();
        }

        return x;
    }

    const Matrix& get_z1() { return m_z1; }
    const Matrix& get_z2() { return m_z2; }
    const Matrix& get_evecs() { return m_evecs; }
};



// [[Rcpp::export]]
List pca_pen_(MapMat S, IntegerVector gr, MapMat x0, double lambda, double gamma,
              double lr, int maxiter = 500, int fan_maxinc = 100, int fan_maxiter = 10,
              double eps = 1e-4, int verbose = 0)
{
    // Dimension of the covariance matrix
    const int n = S.rows();
    const int p = S.cols();
    if(n != p)
        Rcpp::stop("S must be square");

    Matrix linear = Matrix::Constant(p, p, -lambda * gamma * gamma);
    const int ngr = gr.length();
    for(int i = 0; i < ngr; i++)
    {
        const int k = gr[i] - 1;
        linear.col(k) /= gamma;
        linear.row(k) /= gamma;
    }
    linear.noalias() += S;

    PPGOptimizer opt(linear);
    opt.init(x0, fan_maxinc, fan_maxiter, true);

    // Metrics in each iteration
    std::vector<double> err_v;
    std::vector<double> err_z1;
    std::vector<double> err_z2;

    int i = 0;
    for(i = 0; i < maxiter; i++)
    {
        if(verbose > 0)
            Rcpp::Rcout << "iter = " << i << std::endl;

        double diffz1 = opt.update_z1(lr, 0.001 / std::sqrt(i + 1.0), verbose);
        err_z1.push_back(diffz1);

        double diffz2 = opt.update_z2(lr, verbose);
        err_z2.push_back(diffz2);

        bool first_iter = (i == 0);
        double diffev = opt.update_evecs(first_iter);

        if(i > 0)
        {
            err_v.push_back(diffev);

            if(verbose > 0)
                Rcpp::Rcout << "  [info] err_v = " << diffev
                            << ", err_z1 = " << err_z1.back()
                            << ", err_z2 = " << err_z2.back()
                            << std::endl << std::endl;

            if(diffev < eps)
                break;
        }
    }

    const Matrix& x = opt.get_saprse_x(lr);

    return List::create(
        Rcpp::Named("projection") = x,
        Rcpp::Named("evecs") = opt.get_evecs(),
        Rcpp::Named("err_v") = err_v,
        Rcpp::Named("err_z1") = err_z1,
        Rcpp::Named("err_z2") = err_z2,
        Rcpp::Named("niter") = i + 1,
        Rcpp::Named("z1") = opt.get_z1(),
        Rcpp::Named("z2") = opt.get_z2()
    );
}
