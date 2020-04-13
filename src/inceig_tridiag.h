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

    int      m_n;               // Dimension of the matrix
    MatrixXd m_Q;               // Tridiagonal decomposition, x = QTQ'
    VectorXd m_tau;             // Scalar factor for Q
    VectorXd m_diag;            // Diagonal elements of T
    VectorXd m_subdiag;         // Sub-diagonal elements of T

    VectorXd m_dcache;          // Work space
    VectorXd m_lcache;          // Work space
    VectorXd m_ucache;          // Work space

    VectorXd m_evals_lg;        // Largest eigenvalues
    MatrixXd m_evecs_lg;        // Eigenvectors for largest eigenvalues
    VectorXd m_evals_sm;        // Smallest eigenvalues
    MatrixXd m_evecs_sm;        // Eigenvectors for Smallest eigenvalues

    int      m_max_evals_lg;    // Maximum number of largest eigenvalues to be computed
    int      m_max_evals_sm;    // Maximum number of smallest eigenvalues to be computed
    int      m_num_computed_lg; // Number of largest eigenvalues computed
    int      m_num_computed_sm; // Number of smallest eigenvalues computed

    double   m_shift_lg;        // Current shift for the largest eigenvalues
    double   m_shift_sm;        // Current shift for the smallest eigenvalues

    bool     m_mode;            // To compute the largest eigenvalues (true) or the smallest (false)

    // x <- Qx
    //        [  d                  ]
    //        [  e   d              ]
    // Qmat = [  v1  e   d          ]
    //        [  v1  v2  e   d      ]
    //        [  v1  v2  v3  e   d  ]
    // Q = H1 * H2 * ... * Hn-2
    // Hi = I - tau * ui * ui'
    // ui = (0, ..., 0, 1, vi)
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
    // 3. Tridiagonal factorization for the shift-and-invert mode
    inline void init(const RefConstMat& mat, int max_evals_lg, int init_evals_lg, int max_evals_sm, int init_evals_sm)
    {
        // 1. Set the size of the problem
        m_max_evals_lg = max_evals_lg;
        m_num_computed_lg = 0;
        m_shift_lg = 0.0;

        m_max_evals_sm = max_evals_sm;
        m_num_computed_sm = 0;
        m_shift_sm = 0.0;

        m_n = mat.rows();
        m_Q.resize(m_n, m_n);
        m_tau.resize(m_n - 1);
        m_diag.resize(m_n);
        m_subdiag.resize(m_n - 1);

        m_dcache.resize(m_n);
        m_lcache.resize(m_n - 1);
        m_ucache.resize(m_n - 1);

        m_evals_lg.resize(m_max_evals_lg);
        m_evecs_lg.resize(m_n, m_max_evals_lg);
        m_evals_sm.resize(m_max_evals_sm);
        m_evecs_sm.resize(m_n, m_max_evals_sm);

        m_mode = true;

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
        if(init_evals_lg > 0)
        {
            Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, SymTridiag>
            eigs_lg(&op, init_evals_lg, std::min(m_n, std::max(20, init_evals_lg * 2 + 1)));
            eigs_lg.init();
            eigs_lg.compute(1000, 1e-6, Spectra::LARGEST_ALGE);

            // Store computed eigenvalues and eigenvectors
            m_evals_lg.head(init_evals_lg).noalias() = eigs_lg.eigenvalues();
            m_evecs_lg.leftCols(init_evals_lg).noalias() = eigs_lg.eigenvectors();
            m_num_computed_lg += init_evals_lg;
            m_shift_lg = m_evals_lg[m_num_computed_lg - 1] - 1e-6;
        }

        if(init_evals_sm > 0)
        {
            Spectra::SymEigsSolver<double, Spectra::SMALLEST_ALGE, SymTridiag>
            eigs_sm(&op, init_evals_sm, std::min(m_n, std::max(20, init_evals_sm * 2 + 1)));
            eigs_sm.init();
            eigs_sm.compute(1000, 1e-6, Spectra::SMALLEST_ALGE);

            // Store computed eigenvalues and eigenvectors
            m_evals_sm.head(init_evals_sm).noalias() = eigs_sm.eigenvalues();
            m_evecs_sm.leftCols(init_evals_sm).noalias() = eigs_sm.eigenvectors();
            m_num_computed_sm += init_evals_sm;
            m_shift_sm = m_evals_sm[m_num_computed_sm - 1] + 1e-6;
        }
    }

    inline int rows() const { return m_n; }
    inline int cols() const { return m_n; }

    // y = solve(s * I - T, x)
    inline void perform_op(const double* x_in, double* y_out)
    {
        const double shift = m_mode ? m_shift_lg : m_shift_sm;
        // First use the fast solver
        // If not stable (divided-by-zero), use the LAPACK function
        int info = tridiag_shift_solve(
            m_n, m_diag.data(), m_subdiag.data(), x_in, shift, y_out,
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

        m_dcache.array() = shift - m_diag.array();
        m_lcache.noalias() = -m_subdiag;
        m_ucache.noalias() = -m_subdiag;

        int nrhs = 1;
        F77_CALL(dgtsv)(&m_n, &nrhs, m_lcache.data(), m_dcache.data(), m_ucache.data(),
                 y_out, &m_n, &info);
        if(info != 0)
            throw std::runtime_error("failed to compute eigenvalues");
    }

    inline int compute_next_largest(int inc_evals)
    {
        m_mode = true;
        if(m_num_computed_lg + inc_evals > std::min(m_n, m_max_evals_lg))
            throw std::logic_error("maximum number of eigenvalues computed");

        Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, IncrementalEig>
            eigs(this, inc_evals, std::min(m_n, std::max(20, inc_evals * 2 + 1)));
        eigs.init();
        eigs.compute(1000, 1e-6, Spectra::LARGEST_ALGE);

        // Store computed eigenvalues and eigenvectors
        m_evals_lg.segment(m_num_computed_lg, inc_evals).array() = m_shift_lg - 1.0 / eigs.eigenvalues().array();
        m_evecs_lg.block(0, m_num_computed_lg, m_n, inc_evals).noalias() = eigs.eigenvectors();
        m_num_computed_lg += inc_evals;
        m_shift_lg = m_evals_lg[m_num_computed_lg - 1] - 1e-6;

        return eigs.num_operations();
    }

    inline int compute_next_smallest(int inc_evals)
    {
        m_mode = false;
        if(m_num_computed_sm + inc_evals > std::min(m_n, m_max_evals_sm))
            throw std::logic_error("maximum number of eigenvalues computed");

        Spectra::SymEigsSolver<double, Spectra::SMALLEST_ALGE, IncrementalEig>
            eigs(this, inc_evals, std::min(m_n, std::max(20, inc_evals * 2 + 1)));
        eigs.init();
        eigs.compute(1000, 1e-6, Spectra::SMALLEST_ALGE);

        // Store computed eigenvalues and eigenvectors
        m_evals_sm.segment(m_num_computed_sm, inc_evals).array() = m_shift_sm - 1.0 / eigs.eigenvalues().array();
        m_evecs_sm.block(0, m_num_computed_sm, m_n, inc_evals).noalias() = eigs.eigenvectors();
        m_num_computed_sm += inc_evals;
        m_shift_sm = m_evals_sm[m_num_computed_sm - 1] + 1e-6;

        return eigs.num_operations();
    }

    inline void compute_eigenvectors(int num_lg, int num_sm)
    {
        #pragma omp parallel for shared(m_evecs_lg)
        for(int i = 0; i < num_lg; i++)
        {
            apply_Qx(&m_evecs_lg(0, i));
        }
        #pragma omp parallel for shared(m_evecs_sm)
        for(int i = 0; i < num_sm; i++)
        {
            apply_Qx(&m_evecs_sm(0, i));
        }
    }

    const int num_computed_largest() const { return m_num_computed_lg; }
    const int num_computed_smallest() const { return m_num_computed_sm; }
    const VectorXd& largest_eigenvalues() const { return m_evals_lg; }
    const VectorXd& smallest_eigenvalues() const { return m_evals_sm; }
    const MatrixXd& largest_eigenvectors() const { return m_evecs_lg; }
    const MatrixXd& smallest_eigenvectors() const { return m_evecs_sm; }
};


#endif // MARKERPEN_INCEIG_TRIDIAG_H
