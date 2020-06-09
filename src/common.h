#ifndef MARKERPEN_COMMON_H
#define MARKERPEN_COMMON_H

#include <RcppEigen.h>

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using IntVector = Eigen::VectorXi;

using MapMat = Eigen::Map<Matrix>;
using MapVec = Eigen::Map<Vector>;
using MapConstVec = Eigen::Map<const Vector>;
using MapConstMat = Eigen::Map<const Matrix>;

using RefVec = Eigen::Ref<Vector>;
using RefMat = Eigen::Ref<Matrix>;
using RefConstVec = Eigen::Ref<const Vector>;
using RefConstMat = Eigen::Ref<const Matrix>;

using SpMat = Eigen::SparseMatrix<double>;
using MapSpMat = Eigen::Map<SpMat>;
using MapConstSpMat = Eigen::Map<const SpMat>;


#endif  // MARKERPEN_COMMON_H
