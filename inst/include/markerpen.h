#ifndef MARKERPEN_PKG_H
#define MARKERPEN_PKG_H

#include <RcppEigen.h>

typedef Eigen::MatrixXd MatrixXd;
typedef Eigen::VectorXd VectorXd;
typedef Eigen::Map<MatrixXd> MapMat;
typedef Eigen::Map<VectorXd> MapVec;
typedef Eigen::Map<const VectorXd> MapConstVec;
typedef Eigen::Ref<VectorXd> RefVec;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Map<SpMat> MapSpMat;
typedef Eigen::Map<const SpMat> MapConstSpMat;

#endif  // MARKERPEN_PKG_H
