#include "mathbox/decompose.hpp"

namespace math {

template Eigen::MatrixXd LLT<Eigen::MatrixXd>(const Eigen::MatrixBase<Eigen::MatrixXd>&, const LLTDecompositionMethod);
template Eigen::Matrix2d LLT<Eigen::Matrix2d>(const Eigen::MatrixBase<Eigen::Matrix2d>&, const LLTDecompositionMethod);
template Eigen::Matrix3d LLT<Eigen::Matrix3d>(const Eigen::MatrixBase<Eigen::Matrix3d>&, const LLTDecompositionMethod);
template Eigen::Matrix4d LLT<Eigen::Matrix4d>(const Eigen::MatrixBase<Eigen::Matrix4d>&, const LLTDecompositionMethod);

template Eigen::MatrixXd UTU<Eigen::MatrixXd>(const Eigen::MatrixBase<Eigen::MatrixXd>&, const LLTDecompositionMethod);
template Eigen::Matrix2d UTU<Eigen::Matrix2d>(const Eigen::MatrixBase<Eigen::Matrix2d>&, const LLTDecompositionMethod);
template Eigen::Matrix3d UTU<Eigen::Matrix3d>(const Eigen::MatrixBase<Eigen::Matrix3d>&, const LLTDecompositionMethod);
template Eigen::Matrix4d UTU<Eigen::Matrix4d>(const Eigen::MatrixBase<Eigen::Matrix4d>&, const LLTDecompositionMethod);

}
