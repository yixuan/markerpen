#ifndef MARKERPEN_TRIDIAG_H
#define MARKERPEN_TRIDIAG_H

#include <Rcpp.h>

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

class SymTridiag
{
private:
    const int m_n;
    const double* m_diag;
    const double* m_subdiag;

public:
    SymTridiag(int n, const double* diag, const double* subdiag) :
        m_n(n), m_diag(diag), m_subdiag(subdiag)
    {}
    inline int rows() const { return m_n; }
    inline int cols() const { return m_n; }
    inline void perform_op(const double* x_in, double* y_out) const
    {
        tridiag_prod(m_n, m_diag, m_subdiag, x_in, y_out);
    }
};

// Solve (A - s * I) * x = d
// Diagonal:    b[0], ..., b[n-1]
// Subdiagonal: c[0], ..., c[n-2]
// RHS:         d[0], ..., d[n-1]
// Working space - wdc[n-1], wdd[n]. If NULL, they will be created
// Return >0 if there is a divided-by-zero issue
// https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
inline int tridiag_shift_solve(
    int n, const double* b, const double* c, const double* d, double s, double* x,
    double* wdc, double* wdd
)
{
    double* cmod = (wdc == NULL) ? (new double[n - 1]) : wdc;
    double* dmod = (wdd == NULL) ? (new double[n]) : wdd;

    // Tolerance for denominator
    const double eps = 1e-10;

    // Compute new c and d
    double denom = b[0] - s;
    if(std::abs(denom) < eps)
        return 1;
    cmod[0] = c[0] / denom;
    dmod[0] = d[0] / denom;
    for(int i = 1; i < n - 1; i++)
    {
        denom = b[i] - s - c[i - 1] * cmod[i - 1];
        if(std::abs(denom) < eps)
            return 1;
        cmod[i] = c[i] / denom;
        dmod[i] = (d[i] - c[i - 1] * dmod[i - 1]) / denom;
    }
    denom = b[n - 1] - s - c[n - 2] * cmod[n - 2];
    if(std::abs(denom) < eps)
        return 1;
    dmod[n - 1] = (d[n - 1] - c[n - 2] * dmod[n - 2]) / denom;

    // Compute x
    x[n - 1] = dmod[n - 1];
    for(int i = n - 2; i >= 0; i--)
    {
        x[i] = dmod[i] - cmod[i] * x[i + 1];
    }

    if(wdc == NULL)
        delete [] cmod;
    if(wdd == NULL)
        delete [] dmod;

    return 0;
}


#endif // MARKERPEN_TRIDIAG_H
