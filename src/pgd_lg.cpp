#include <RcppEigen.h>
#include "common.h"

using Rcpp::List;

inline void proj_pos_simplex(MapVec x, double a)
{
    const int n = x.size();

    // Sort x in decreasing order
    VectorXd u = x;
    std::sort(u.data(), u.data() + n, std::greater<double>());

    VectorXd cumu(n);
    cumu[0] = u[0];
    for(int i = 1; i < n; i++)
        cumu[i] = cumu[i - 1] + u[i];

    int j = n;
    double lambda = 0.0;
    for(; j > 0; j--)
    {
        lambda = (a - cumu[j - 1]) / j;
        if(u[j - 1] + lambda > 0.0)
            break;
    }

    x.array() = (x.array() + lambda).max(0.0);
}

// [[Rcpp::export]]
List pgd_lg(MapMat X, MapMat W, MapMat F0, MapVec S0, MapMat G0,
            double lambda, double lr, int niter,
            bool upd_F = true, bool upd_S = true, bool upd_G = true,
            double log_eps = 1e-6,
            int verbose_freq = 10)
{
    const int n = X.rows();
    const int p = X.cols();
    const int nc = S0.size();

    MatrixXd Y = (X.array() + log_eps).log().matrix();
    MatrixXd F = F0;
    VectorXd S = S0;
    MatrixXd G = G0;

    std::vector<double> loss;
    MatrixXd GS(p, nc), R(n, p), logR(n, p), dR(n, p), dF(n, nc), dG(p, nc), A(n * p, nc);
    VectorXd dS(nc), Frow(nc);

    for(int i = 0; i < niter; i++)
    {
        GS.noalias() = G * S.asDiagonal();
        R.noalias() = F * GS.transpose();
        logR.array() = (R.array() + log_eps).log();

        double obj1 = (Y - logR).squaredNorm();
        double obj2 = W.cwiseProduct(G.cwiseAbs2()).sum();
        double obj = obj1 + lambda * obj2;
        loss.push_back(obj);
        if(i % verbose_freq == 0)
            ::Rprintf("k = %d, obj1 = %f, obj2 = %f, obj = %f\n", i, obj1, obj2, obj);

        dR.array() = 2.0 * (logR.array() - Y.array()) / (R.array() + log_eps);
        dF.noalias() = dR * GS;

        for(int j = 0; j < nc; j++)
        {
            MapMat(A.col(j).data(), n, p).noalias() = F.col(j) * G.col(j).transpose();
        }
        dS.noalias() = A.transpose() * MapVec(dR.data(), n * p);

        dG.noalias() = dR.transpose() * F * S.asDiagonal() + 2.0 * lambda * G.cwiseProduct(W);
        if(i % verbose_freq == 0)
            ::Rprintf("===> ||dF|| = %f, ||dS|| = %f, ||dG|| = %f\n\n",
                      dF.norm(), dS.norm(), dG.norm());

        if(upd_F)
        {
            F.noalias() -= lr * dF;
            for(int j = 0; j < n; j++)
            {
                Frow.noalias() = F.row(j).transpose();
                proj_pos_simplex(MapVec(Frow.data(), nc), 1.0);
                F.row(j).noalias() = Frow;
            }
        }

        if(upd_S)
        {
            S.noalias() -= lr * dS;
            S.array() = S.array().max(0.0);
        }

        if(upd_G)
        {
            G.noalias() -= lr * dG;
            for(int j = 0; j < nc; j++)
            {
                proj_pos_simplex(MapVec(G.col(j).data(), p), 1.0);
            }
        }
    }

    return List::create(
        Rcpp::Named("loss") = loss,
        Rcpp::Named("F") = F,
        Rcpp::Named("S") = S,
        Rcpp::Named("G") = G
    );
}

