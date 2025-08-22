#include "mathbox/geometry.hpp"

namespace math {

template Eigen::Isometry3d change_relative_transform_frame<double>(const Eigen::Isometry3d&, const Eigen::Isometry3d&);

template Eigen::Isometry3d change_relative_transform_frame<double>(const Eigen::Isometry3d&, const Eigen::Isometry3d&,
        const Eigen::Isometry3d&);

template Eigen::Matrix<double, 6, 6> change_tf_covariance_frame<double>(const Eigen::Matrix<double, 6, 6>&,
        const Eigen::Isometry3d&, const bool);

template Eigen::Matrix<double, 6, 1> change_twist_reference_frame<double>(const Eigen::Isometry3d&,
        const Eigen::Matrix<double, 6, 1>&, const bool);

template Eigen::Matrix<double, 3, 3> compose_transform_covariance<double, 2>(const Eigen::Matrix<double, 3, 3>&,
        const Eigen::Matrix<double, 3, 3>&, const Eigen::Isometry2d&, const Eigen::Matrix<double, 3, 3>&, const bool);
template Eigen::Matrix<double, 6, 6> compose_transform_covariance<double, 3>(const Eigen::Matrix<double, 6, 6>&,
        const Eigen::Matrix<double, 6, 6>&, const Eigen::Isometry3d&, const Eigen::Matrix<double, 6, 6>&, const bool);

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

template Eigen::Matrix<double, 3, 3> transform_adjoint<double, 2>(const Eigen::Isometry2d& transform, const bool);
template Eigen::Matrix<double, 6, 6> transform_adjoint<double, 3>(const Eigen::Isometry3d& transform, const bool);

}
