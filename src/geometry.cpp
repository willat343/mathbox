#include "mathbox/geometry.hpp"

namespace math {

template Eigen::Isometry3d change_relative_transform_frame<double>(const Eigen::Isometry3d&, const Eigen::Isometry3d&);

template Eigen::Isometry3d change_relative_transform_frame<double>(const Eigen::Isometry3d&, const Eigen::Isometry3d&,
        const Eigen::Isometry3d&);

template Eigen::Matrix<double, 6, 6> change_tf_covariance_frame_tr<double>(const Eigen::Matrix<double, 6, 6>&,
        const Eigen::Isometry3d&);

template Eigen::Matrix<double, 6, 1> change_twist_reference_frame_tr<double>(const Eigen::Isometry3d&,
        const Eigen::Matrix<double, 6, 1>&);

template Eigen::Matrix<double, 3, 3> compose_transform_covariance_tr<double, 2>(const Eigen::Matrix<double, 3, 3>&,
        const Eigen::Matrix<double, 3, 3>&, const Eigen::Isometry2d&, const Eigen::Matrix<double, 3, 3>&);
template Eigen::Matrix<double, 6, 6> compose_transform_covariance_tr<double, 3>(const Eigen::Matrix<double, 6, 6>&,
        const Eigen::Matrix<double, 6, 6>&, const Eigen::Isometry3d&, const Eigen::Matrix<double, 6, 6>&);

template Eigen::Matrix<double, 6, 1> compute_constant_rates<double>(const typename Eigen::Isometry3d& pose_1,
        const typename Eigen::Isometry3d&, const double);

template Eigen::Isometry3d glerp<double>(const Eigen::Isometry3d& T_0, const Eigen::Isometry3d& T_1, const double);

template Eigen::Isometry2d relative_transform<double, 2>(const typename Eigen::Isometry2d&,
        const typename Eigen::Isometry2d&);
template Eigen::Isometry3d relative_transform<double, 3>(const typename Eigen::Isometry3d&,
        const typename Eigen::Isometry3d&);

template Eigen::Matrix2d rotate_point_covariance<Eigen::Matrix2d>(const Eigen::MatrixBase<Eigen::Matrix2d>&,
        const Eigen::MatrixBase<Eigen::Matrix2d>&);
template Eigen::Matrix3d rotate_point_covariance<Eigen::Matrix3d>(const Eigen::MatrixBase<Eigen::Matrix3d>&,
        const Eigen::MatrixBase<Eigen::Matrix3d>&);

#ifndef DOXYGEN_EXCLUDE
template Eigen::Vector2d so_cross<Eigen::Vector<double, 1>, Eigen::Vector2d>(
        const Eigen::MatrixBase<Eigen::Vector<double, 1>>&, const Eigen::MatrixBase<Eigen::Vector2d>&);
template Eigen::Vector3d so_cross<Eigen::Vector3d, Eigen::Vector3d>(const Eigen::MatrixBase<Eigen::Vector3d>&,
        const Eigen::MatrixBase<Eigen::Vector3d>&);
#endif

template Eigen::Vector<double, 1> so_from_skew<Eigen::Matrix2d>(const Eigen::MatrixBase<Eigen::Matrix2d>&);
template Eigen::Vector3d so_from_skew<Eigen::Matrix3d>(const Eigen::MatrixBase<Eigen::Matrix3d>&);

template Eigen::Matrix2d so_skew<Eigen::Vector<double, 1>>(const Eigen::MatrixBase<Eigen::Vector<double, 1>>&);
template Eigen::Matrix3d so_skew<Eigen::Vector3d>(const Eigen::MatrixBase<Eigen::Vector3d>&);

template Eigen::Matrix<double, 3, 3> adjoint_SE_tr<double, 2>(const Eigen::Isometry2d& transform);
template Eigen::Matrix<double, 6, 6> adjoint_SE_tr<double, 3>(const Eigen::Isometry3d& transform);

}
