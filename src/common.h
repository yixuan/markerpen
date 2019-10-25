#ifndef MARKERPEN_COMMON_H
#define MARKERPEN_COMMON_H

#include <RcppEigen.h>

typedef Eigen::MatrixXd MatrixXd;
typedef Eigen::VectorXd VectorXd;

typedef Eigen::Map<MatrixXd> MapMat;
typedef Eigen::Map<VectorXd> MapVec;
typedef Eigen::Map<const VectorXd> MapConstVec;
typedef Eigen::Map<const MatrixXd> MapConstMat;

typedef Eigen::Ref<VectorXd> RefVec;
typedef Eigen::Ref<MatrixXd> RefMat;
typedef Eigen::Ref<const VectorXd> RefConstVec;
typedef Eigen::Ref<const MatrixXd> RefConstMat;

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Map<SpMat> MapSpMat;
typedef Eigen::Map<const SpMat> MapConstSpMat;


#endif  // MARKERPEN_COMMON_H
