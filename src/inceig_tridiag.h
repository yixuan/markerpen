#ifndef MARKERPEN_INCEIG_TRIDIAG_H
#define MARKERPEN_INCEIG_TRIDIAG_H

#include "common.h"
#include "tridiag.h"
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include "walltime.h"

#ifdef EIGEN_USE_BLAS
#define F77_CALL(x)	x ## _
#define F77_NAME(x) F77_CALL(x)
#define La_extern extern
extern "C" {

La_extern void
F77_NAME(dsytrd)(const char* uplo, const int* n,
                 double* a, const int* lda,
                 double* d, double* e, double* tau,
                 double* work, const int* lwork, int* info);

La_extern void
F77_NAME(dgtsv)(const int* n, const int* nrhs,
                double* dl, double* d, double* du,
                double* b, const int* ldb, int* info);

}
#else
#include <R_ext/Lapack.h>
#endif


class IncrementalEig
{
private:
    typedef Eigen::Map<VectorXd> MapVec;
    typedef Eigen::Map<const VectorXd> MapConstVec;
    typedef Eigen::Map<const MatrixXd> MapConstMat;
    typedef Eigen::Ref<const MatrixXd> RefConstMat;

    int           m_n;
    MatrixXd      m_Q;
    VectorXd      m_tau;
    VectorXd      m_diag;
    VectorXd      m_subdiag;

    VectorXd      m_dcache;
    VectorXd      m_lcache;
    VectorXd      m_ucache;

    VectorXd      m_evals;
    MatrixXd      m_evecs;

    int           m_max_evals;
    int           m_init_evals;
    int           m_num_computed;

    double        m_shift;

    inline void apply_Qx(double* xptr) const
    {

        int vlen = 1;
        for(int i = m_n - 3; i >= 0; i--, vlen++)
        {
            MapVec xtail(xptr + i + 2, vlen);
            MapConstVec v(&m_Q(i + 2, i), vlen);
            const double vx = v.dot(xtail) + xptr[i + 1];
            const double scale = vx * m_tau[i];
            xtail.noalias() -= scale * v;
            xptr[i + 1] -= scale;
        }
    }

public:
    IncrementalEig() {}

    // 1. Set the size of the problem
    // 2. Compute initial `init_evals` eigenvalues
    // 3. Cholesky factorization for the shift-and-invert mode
    inline void init(const RefConstMat& mat, int max_evals, int init_evals)
    {
        // 1. Set the size of the problem
        m_max_evals = max_evals;
        m_init_evals = init_evals;
        m_num_computed = 0;
        m_shift = 0.0;

        m_n = mat.rows();
        m_Q.resize(m_n, m_n);
        m_tau.resize(m_n - 1);
        m_diag.resize(m_n);
        m_subdiag.resize(m_n - 1);

        m_dcache.resize(m_n);
        m_lcache.resize(m_n - 1);
        m_ucache.resize(m_n - 1);

        m_evals.resize(m_max_evals);
        m_evecs.resize(m_n, m_max_evals);

        // 2. Tridiagonalization
        char uplo = 'L';
        int lwork = -1, info;
        double blocksize;

        m_Q.noalias() = mat;
        F77_CALL(dsytrd)(&uplo, &m_n,
                 m_Q.data(), &m_n,
                 m_diag.data(), m_subdiag.data(), m_tau.data(),
                 &blocksize, &lwork, &info);

        lwork = int(blocksize);
        std::vector<double> work;
        work.resize(lwork);

        F77_CALL(dsytrd)(&uplo, &m_n,
                 m_Q.data(), &m_n,
                 m_diag.data(), m_subdiag.data(), m_tau.data(),
                 &work[0], &lwork, &info);

        // 3. Compute initial `init_evals` eigenvalues
        SymTridiag op(m_n, m_diag.data(), m_subdiag.data());
        Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, SymTridiag>
            eigs(&op, init_evals, std::min(m_n, std::max(20, init_evals * 2 + 1)));
        eigs.init();
        eigs.compute(1000, 1e-6);

        // Store computed eigenvalues and eigenvectors
        m_evals.head(init_evals).noalias() = eigs.eigenvalues();
        m_evecs.leftCols(init_evals).noalias() = eigs.eigenvectors();
        m_num_computed += init_evals;
        m_shift = m_evals[m_num_computed - 1] - 1e-6;
    }

    inline int rows() const { return m_n; }
    inline int cols() const { return m_n; }

    // y = solve(s * I - T, x)
    inline void perform_op(const double* x_in, double* y_out)
    {
        // First use the fast solver
        // If not stable (divided-by-zero), use the LAPACK function
        int info = tridiag_shift_solve(
            m_n, m_diag.data(), m_subdiag.data(), x_in, m_shift, y_out,
            m_lcache.data(), m_dcache.data()
        );
        if(info == 0)
        {
            // Negate y
            std::transform(y_out, y_out + m_n, y_out, std::negate<double>());
            return;
        }

        // Otherwise, use LAPACK
        std::copy(x_in, x_in + m_n, y_out);

        m_dcache.array() = m_shift - m_diag.array();
        m_lcache.noalias() = -m_subdiag;
        m_ucache.noalias() = -m_subdiag;

        int nrhs = 1;
        F77_CALL(dgtsv)(&m_n, &nrhs, m_lcache.data(), m_dcache.data(), m_ucache.data(),
                 y_out, &m_n, &info);
        if(info != 0)
            throw std::runtime_error("failed to compute eigenvalues");
    }

    inline int compute_next(int inc_evals)
    {
        if(m_num_computed + inc_evals > m_max_evals)
            throw std::logic_error("maximum number of eigenvalues computed");

        Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, IncrementalEig>
            eigs(this, inc_evals, std::min(m_n, std::max(20, inc_evals * 2 + 1)));
        eigs.init();
        eigs.compute(1000, 1e-6);

        // Store computed eigenvalues and eigenvectors
        m_evals.segment(m_num_computed, inc_evals).array() = m_shift - 1.0 / eigs.eigenvalues().array();
        m_evecs.block(0, m_num_computed, m_n, inc_evals).noalias() = eigs.eigenvectors();
        m_num_computed += inc_evals;
        m_shift = m_evals[m_num_computed - 1] - 1e-6;

        return eigs.num_operations();
    }

    inline void compute_eigenvectors()
    {
        #pragma omp parallel for shared(m_num_computed, m_evecs)
        for(int i = 0; i < m_num_computed; i++)
        {
            apply_Qx(&m_evecs(0, i));
        }
    }

    const int num_computed() const { return m_num_computed; }
    const VectorXd& eigenvalues() const { return m_evals; }
    const MatrixXd& eigenvectors() const { return m_evecs; }
};


#endif // MARKERPEN_INCEIG_TRIDIAG_H
