#ifndef MARKERPEN_INCEIG_TRIDIAG_LAPACK_H
#define MARKERPEN_INCEIG_TRIDIAG_LAPACK_H

#include "common.h"
#include "walltime.h"

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

La_extern void
F77_NAME(dstegr)(const char* jobz, const char* range, const int* n, double* d,
                 double* e, const double* vl, const double* vu,
                 const int* il, const int* iu, const double* abstol,
                 int* m, double* w, double* z, const int* ldz, int* isuppz,
                 double* work, const int* lwork, int* iwork, const int* liwork, int* info);

}


class IncrementalEig
{
private:
    const int  m_n;               // Dimension of the matrix
    Matrix     m_Q;               // Tridiagonal decomposition, x = QTQ'
    Vector     m_tau;             // Scalar factor for Q
    Vector     m_diag;            // Diagonal elements of T
    Vector     m_subdiag;         // Sub-diagonal elements of T

    Vector     m_evals_lg;        // Largest eigenvalues
    Matrix     m_evecs_lg;        // Eigenvectors for largest eigenvalues
    Vector     m_evals_sm;        // Smallest eigenvalues
    Matrix     m_evecs_sm;        // Eigenvectors for Smallest eigenvalues

    int        m_max_evals_lg;    // Maximum number of largest eigenvalues to be computed
    int        m_max_evals_sm;    // Maximum number of smallest eigenvalues to be computed
    int        m_num_computed_lg; // Number of largest eigenvalues computed
    int        m_num_computed_sm; // Number of smallest eigenvalues computed

    Vector     m_workd;
    Vector     m_worke;
    Vector     m_workevals;
    IntVector  m_suppz;
    std::vector<double> m_workarr;
    std::vector<int>    m_iworkarr;

    inline void compute_eigenpair(int il, int iu, double* evals, double* evecs, bool sort_desc = true)
    {
        constexpr char jobz = 'V';
        constexpr char range = 'I';
        int n = m_n;
        m_workd.noalias() = m_diag;
        m_worke.head(m_n - 1).noalias() = m_subdiag;
        int nevals, lwork = -1, liwork = -1;
        double blocksize;
        int iblocksize;
        int info;
        F77_CALL(dstegr)(&jobz, &range, &n, m_workd.data(), m_worke.data(), NULL, NULL,
                 &il, &iu, NULL,
                 &nevals, m_workevals.data(), evecs, &n, m_suppz.data(),
                 &blocksize, &lwork, &iblocksize, &liwork, &info);
        if(info != 0)
            throw std::runtime_error("computing eigenpair failed");

        lwork = int(blocksize);
        liwork = iblocksize;

        m_workarr.reserve(lwork);
        m_iworkarr.reserve(liwork);

        F77_CALL(dstegr)(&jobz, &range, &n, m_workd.data(), m_worke.data(), NULL, NULL,
                 &il, &iu, NULL,
                 &nevals, m_workevals.data(), evecs, &n, m_suppz.data(),
                 &m_workarr[0], &lwork, &m_iworkarr[0], &liwork, &info);
        if(info != 0)
            throw std::runtime_error("computing eigenpair failed");

        if(sort_desc)
        {
            for(int i = 0; i < nevals; i++)
                evals[i] = m_workevals[nevals - i - 1];
            for(int i = 0; i < nevals / 2; i++)
            {
                std::swap_ranges(evecs + i * m_n, evecs + (i + 1) * m_n,
                                 evecs + (nevals - i - 1) * m_n);
            }
        } else {
            std::copy(m_workevals.data(), m_workevals.data() + nevals, evals);
        }
    }

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
        m_workd(m_n), m_worke(m_n), m_workevals(m_n), m_suppz(2 * m_n)
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

        m_max_evals_sm = max_evals_sm;
        m_num_computed_sm = 0;

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
        if(init_evals_lg > 0)
        {
            compute_eigenpair(m_n - init_evals_lg + 1, m_n,
                              m_evals_lg.data(), m_evecs_lg.data(), true);
            m_num_computed_lg += init_evals_lg;
        }
        if(init_evals_sm > 0)
        {
            compute_eigenpair(1, init_evals_sm,
                              m_evals_sm.data(), m_evecs_sm.data(), false);
            m_num_computed_sm += init_evals_sm;
        }
    }

    inline int compute_next_largest(int inc_evals)
    {
        if(m_num_computed_lg + inc_evals > std::min(m_n, m_max_evals_lg))
            throw std::logic_error("maximum number of eigenvalues computed");

        compute_eigenpair(m_n - m_num_computed_lg - inc_evals + 1, m_n - m_num_computed_lg,
                          &m_evals_lg[m_num_computed_lg], &m_evecs_lg(0, m_num_computed_lg), true);
        m_num_computed_lg += inc_evals;

        return 0;
    }

    inline int compute_next_smallest(int inc_evals)
    {
        if(m_num_computed_sm + inc_evals > std::min(m_n, m_max_evals_sm))
            throw std::logic_error("maximum number of eigenvalues computed");

        compute_eigenpair(m_num_computed_sm + 1, m_num_computed_sm + inc_evals,
                          &m_evals_sm[m_num_computed_sm], &m_evecs_sm(0, m_num_computed_sm), false);
        m_num_computed_sm += inc_evals;

        return 0;
    }

    inline void compute_eigenvectors(int num_lg, int num_sm) noexcept
    {
#if defined(_OPENMP)
        #pragma omp parallel for
#endif
        for(int i = 0; i < num_lg; i++)
        {
            apply_Qx(&m_evecs_lg(0, i));
        }
#if defined(_OPENMP)
        #pragma omp parallel for
#endif
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


#endif // MARKERPEN_INCEIG_TRIDIAG_LAPACK_H
