#include "mathbox/matrix_operations.hpp"

namespace math {

template Eigen::Matrix3d skew_symmetric_cross(const Eigen::MatrixBase<Eigen::Vector3d>& v);

template Eigen::Matrix<double, 3, 2> orthonormal_complement(const Eigen::MatrixBase<Eigen::Vector3d>& unit_vector);

}
