#ifndef MARKERPEN_INCEIG_TRIDIAG_SPECTRA_H
#define MARKERPEN_INCEIG_TRIDIAG_SPECTRA_H

#include "common.h"
#include "tridiag.h"
#include "walltime.h"
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>

#ifndef F77_CALL
#define F77_CALL(x)	x ## _
#endif

#ifndef F77_NAME
#define F77_NAME(x) F77_CALL(x)
#endif

#ifndef La_extern
#define La_extern extern
#endif

extern "C" {

La_extern void
F77_NAME(dsytrd)(const char* uplo, const int* n,
                 double* a, const int* lda,
                 double* d, double* e, double* tau,
                 double* work, const int* lwork, int* info);

}


class IncrementalEig
{
private:
    const int  m_n;               // Dimension of the matrix
    Matrix     m_Q;               // Tridiagonal decomposition, x = QTQ'
    Vector     m_tau;             // Scalar factor for Q
    Vector     m_diag;            // Diagonal elements of T
    Vector     m_subdiag;         // Sub-diagonal elements of T
    SymTridiag m_tridiagop;       // Operator for the tridiagonal matrix

    Vector     m_evals_lg;        // Largest eigenvalues
    Matrix     m_evecs_lg;        // Eigenvectors for largest eigenvalues
    Vector     m_evals_sm;        // Smallest eigenvalues
    Matrix     m_evecs_sm;        // Eigenvectors for Smallest eigenvalues

    int        m_max_evals_lg;    // Maximum number of largest eigenvalues to be computed
    int        m_max_evals_sm;    // Maximum number of smallest eigenvalues to be computed
    int        m_num_computed_lg; // Number of largest eigenvalues computed
    int        m_num_computed_sm; // Number of smallest eigenvalues computed

    double     m_shift_lg;        // Current shift for the largest eigenvalues
    double     m_shift_sm;        // Current shift for the smallest eigenvalues

    // x <- Qx
    //        [  d                  ]
    //        [  e   d              ]
    // Qmat = [  v1  e   d          ]
    //        [  v1  v2  e   d      ]
    //        [  v1  v2  v3  e   d  ]
    // Q = H1 * H2 * ... * Hn-2
    // Hi = I - tau * ui * ui'
    // ui = (0, ..., 0, 1, vi)
    inline void apply_Qx(double* xptr) const noexcept
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
    IncrementalEig(int n):
        m_n(n), m_Q(n, n), m_tau(n - 1), m_diag(n), m_subdiag(n - 1),
        m_tridiagop(n, m_diag.data(), m_subdiag.data())
    {}

    // 1. Set the size of the problem
    // 2. Compute initial `init_evals` eigenvalues
    // 3. Tridiagonal factorization for the shift-and-invert mode
    inline void init(const RefConstMat& mat, int max_evals_lg, int init_evals_lg, int max_evals_sm, int init_evals_sm)
    {
        if(mat.rows() != m_n || mat.cols() != m_n)
            throw std::invalid_argument("matrix dimensions do not match");

        // 1. Set the size of the problem
        m_max_evals_lg = max_evals_lg;
        m_num_computed_lg = 0;
        m_shift_lg = 0.0;

        m_max_evals_sm = max_evals_sm;
        m_num_computed_sm = 0;
        m_shift_sm = 0.0;

        m_evals_lg.resize(m_max_evals_lg);
        m_evecs_lg.resize(m_n, m_max_evals_lg);
        m_evals_sm.resize(m_max_evals_sm);
        m_evecs_sm.resize(m_n, m_max_evals_sm);

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
        work.reserve(lwork);

        F77_CALL(dsytrd)(&uplo, &m_n,
                 m_Q.data(), &m_n,
                 m_diag.data(), m_subdiag.data(), m_tau.data(),
                 &work[0], &lwork, &info);

        // 3. Compute initial `init_evals` eigenvalues
        m_tridiagop.set_mode(OpMode::Prod);
        if(init_evals_lg > 0)
        {
            Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, SymTridiag>
            eigs_lg(&m_tridiagop, init_evals_lg, std::min(m_n, std::max(20, init_evals_lg * 2 + 1)));
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
            eigs_sm(&m_tridiagop, init_evals_sm, std::min(m_n, std::max(20, init_evals_sm * 2 + 1)));
            eigs_sm.init();
            eigs_sm.compute(1000, 1e-6, Spectra::SMALLEST_ALGE);

            // Store computed eigenvalues and eigenvectors
            m_evals_sm.head(init_evals_sm).noalias() = eigs_sm.eigenvalues();
            m_evecs_sm.leftCols(init_evals_sm).noalias() = eigs_sm.eigenvectors();
            m_num_computed_sm += init_evals_sm;
            m_shift_sm = m_evals_sm[m_num_computed_sm - 1] + 1e-6;
        }
    }

    inline int compute_next_largest(int inc_evals)
    {
        if(m_num_computed_lg + inc_evals > std::min(m_n, m_max_evals_lg))
            throw std::logic_error("maximum number of eigenvalues computed");

        m_tridiagop.set_mode(OpMode::ShiftSolve);
        bool success = m_tridiagop.factorize(m_shift_lg);
        if(!success)
            throw std::runtime_error("failed to compute eigenvalues");

        Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, SymTridiag>
            eigs(&m_tridiagop, inc_evals, std::min(m_n, std::max(20, inc_evals * 2 + 1)));
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
        if(m_num_computed_sm + inc_evals > std::min(m_n, m_max_evals_sm))
            throw std::logic_error("maximum number of eigenvalues computed");

        m_tridiagop.set_mode(OpMode::ShiftSolve);
        bool success = m_tridiagop.factorize(m_shift_sm);
        if(!success)
            throw std::runtime_error("failed to compute eigenvalues");

        Spectra::SymEigsSolver<double, Spectra::SMALLEST_ALGE, SymTridiag>
            eigs(&m_tridiagop, inc_evals, std::min(m_n, std::max(20, inc_evals * 2 + 1)));
        eigs.init();
        eigs.compute(1000, 1e-6, Spectra::SMALLEST_ALGE);

        // Store computed eigenvalues and eigenvectors
        m_evals_sm.segment(m_num_computed_sm, inc_evals).array() = m_shift_sm - 1.0 / eigs.eigenvalues().array();
        m_evecs_sm.block(0, m_num_computed_sm, m_n, inc_evals).noalias() = eigs.eigenvectors();
        m_num_computed_sm += inc_evals;
        m_shift_sm = m_evals_sm[m_num_computed_sm - 1] + 1e-6;

        return eigs.num_operations();
    }

    inline void compute_eigenvectors(int num_lg, int num_sm) noexcept
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
    const Vector& largest_eigenvalues() const { return m_evals_lg; }
    const Vector& smallest_eigenvalues() const { return m_evals_sm; }
    const Matrix& largest_eigenvectors() const { return m_evecs_lg; }
    const Matrix& smallest_eigenvectors() const { return m_evecs_sm; }
};


#endif // MARKERPEN_INCEIG_TRIDIAG_SPECTRA_H
