#ifndef MARKERPEN_TRIDIAG_H
#define MARKERPEN_TRIDIAG_H

#include <cmath>
#include <cstdlib>
#include <vector>

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
F77_NAME(dgttrf)(const int* n, double* dl, double* d,
                 double* du, double* du2, int* ipiv, int* info);

La_extern void
F77_NAME(dgttrs)(const char* trans, const int* n, const int* nrhs,
                 double* dl, double* d, double* du, double* du2,
                 int* ipiv, double* b, const int* ldb, int* info);

}

// y = A * x
// Diagonal:    b[0], ..., b[n-1]
// Subdiagonal: c[0], ..., c[n-2]
inline void tridiag_prod(int n, const double* b, const double* c, const double* x, double* y)
{
    y[0] = b[0] * x[0] + c[0] * x[1];
    y[n - 1] = c[n - 2] * x[n - 2] + b[n - 1] * x[n - 1];
    for(int i = 1; i < n - 1; i++)
    {
        y[i] = c[i - 1] * x[i - 1] + b[i] * x[i] + c[i] * x[i + 1];
    }
}

// Factorize A to solve (A - s * I) * x = d
// Diagonal:    b[0], ..., b[n-1]
// Subdiagonal: c[0], ..., c[n-2]
// Working space cmod[n-1], denom[n]
// Return true if successful, otherwise there is a divided-by-zero issue
// https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
inline bool tridiag_fac(int n, const double* b, const double* c, double s,
                        double* cmod, double* denom)
{
    // Tolerance for denominator
    constexpr double eps = 1e-10;

    // Compute the modified c and denominators
    double dtest = b[0] - s;
    if(std::abs(dtest) < eps)
        return false;
    denom[0] = 1.0 / dtest;
    cmod[0] = c[0] * denom[0];
    for(int i = 1; i < n - 1; i++)
    {
        dtest = b[i] - s - c[i - 1] * cmod[i - 1];
        if(std::abs(dtest) < eps)
            return false;
        denom[i] = 1.0 / dtest;
        cmod[i] = c[i] * denom[i];
    }
    dtest = b[n - 1] - s - c[n - 2] * cmod[n - 2];
    denom[n - 1] = 1.0 / dtest;
    return std::abs(dtest) >= eps;
}

// Solve (A - s * I) * x = d based on the factorization result
inline void tridiag_solve(int n, const double* c, const double* cmod, const double* denom,
                          const double* d, double* x)
{
    // Use x to hold the modified values of d
    x[0] = d[0] * denom[0];
    for(int i = 1; i < n; i++)
    {
        x[i] = (d[i] - c[i - 1] * x[i - 1]) * denom[i];
    }
    // Back substitution
    for(int i = n - 2; i >= 0; i--)
    {
        x[i] -= cmod[i] * x[i + 1];
    }
}

enum class OpMode { Prod, ShiftSolve };
enum class ShiftSolver { Fast, Lapack };

class SymTridiag
{
public:
    using Scalar = double;

private:
    const int m_n;
    const double* m_diag;
    const double* m_subdiag;

    // Operation mode
    OpMode m_mode;
    // Solver
    ShiftSolver m_solver;

    // Working spaces
    mutable std::vector<double> m_dcache;
    mutable std::vector<double> m_lcache;
    mutable std::vector<double> m_ucache;
    mutable std::vector<double> m_u2cache;
    mutable std::vector<int> m_icache;

public:
    SymTridiag(int n, const double* diag, const double* subdiag) :
        m_n(n), m_diag(diag), m_subdiag(subdiag),
        m_mode(OpMode::Prod), m_solver(ShiftSolver::Fast)
    {
        // By default we don't allocate memory for caches at this moment,
        // since they are only needed when factorize() is called
    }

    inline int rows() const { return m_n; }
    inline int cols() const { return m_n; }
    inline void set_mode(OpMode mode) noexcept { m_mode = mode; }

    inline bool factorize(double shift) noexcept
    {
        // First use the fast solver
        // If not stable (divided-by-zero), use the LAPACK function
        //
        // Here we need to allocate space for l and d caches
        m_dcache.reserve(m_n);
        m_lcache.reserve(m_n - 1);
        bool success = tridiag_fac(m_n, m_diag, m_subdiag, shift,
                                   &m_lcache[0], &m_dcache[0]);
        if(success)
        {
            m_solver = ShiftSolver::Fast;
            return true;
        }

        // If not successful, use LU decomposition
        //
        // Allocate memory for caches
        m_ucache.reserve(m_n - 1);
        m_u2cache.reserve(m_n - 2);
        m_icache.reserve(m_n);
        // Make copies of coefficients
        std::transform(m_diag, m_diag + m_n, m_dcache.begin(),
                       [shift](double x){ return shift - x; });
        std::transform(m_subdiag, m_subdiag + m_n - 1, m_lcache.begin(),
                       std::negate<double>());
        std::copy(m_lcache.begin(), m_lcache.end(), m_ucache.begin());

        int info;
        F77_CALL(dgttrf)(&m_n, &m_lcache[0], &m_dcache[0], &m_ucache[0],
                         &m_u2cache[0], &m_icache[0], &info);
        m_solver = ShiftSolver::Lapack;
        return info == 0;
    }

    inline void perform_op(const double* x_in, double* y_out) const noexcept
    {
        // Computing product
        if(m_mode == OpMode::Prod)
        {
            tridiag_prod(m_n, m_diag, m_subdiag, x_in, y_out);
            return;
        }

        // Fast shift solver
        if(m_solver == ShiftSolver::Fast)
        {
            tridiag_solve(m_n, m_subdiag, &m_lcache[0], &m_dcache[0], x_in, y_out);
            // Negate y
            std::transform(y_out, y_out + m_n, y_out, std::negate<double>());
            return;
        }

        // LU-based shift solver
        constexpr char trans = 'N';
        constexpr int nrhs = 1;
        int info;
        std::copy(x_in, x_in + m_n, y_out);
        F77_CALL(dgttrs)(&trans, &m_n, &nrhs,
                         &m_lcache[0], &m_dcache[0], &m_ucache[0], &m_u2cache[0],
                         &m_icache[0], y_out, &m_n, &info);
    }
};


#endif // MARKERPEN_TRIDIAG_H
