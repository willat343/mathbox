#ifndef MATHBOX_IMPL_GEOMETRY_HPP
#define MATHBOX_IMPL_GEOMETRY_HPP

#include "mathbox/geometry.hpp"
#include "mathbox/matrix_operations.hpp"

namespace math {

template<typename Scalar>
inline Eigen::Transform<Scalar, 3, Eigen::Isometry> change_relative_transform_frame(
        const typename Eigen::Transform<Scalar, 3, Eigen::Isometry>& relative_transform_A,
        const typename Eigen::Transform<Scalar, 3, Eigen::Isometry>& rigid_transform_B_A) {
    return change_relative_transform_frame(relative_transform_A, rigid_transform_B_A, rigid_transform_B_A.inverse());
}

template<typename Scalar>
inline Eigen::Transform<Scalar, 3, Eigen::Isometry> change_relative_transform_frame(
        const typename Eigen::Transform<Scalar, 3, Eigen::Isometry>& relative_transform_A,
        const typename Eigen::Transform<Scalar, 3, Eigen::Isometry>& rigid_transform_B_A,
        const typename Eigen::Transform<Scalar, 3, Eigen::Isometry>& rigid_transform_A_B) {
    return rigid_transform_B_A * relative_transform_A * rigid_transform_A_B;
}

template<typename Scalar>
inline Eigen::Matrix<Scalar, 6, 6> change_tf_covariance_frame(const typename Eigen::Matrix<Scalar, 6, 6>& covariance_A,
        const typename Eigen::Matrix<Scalar, 6, 6>& transform_adjoint_B_A) {
    return transform_adjoint_B_A * covariance_A * transform_adjoint_B_A.transpose();
}

template<typename Scalar>
inline Eigen::Matrix<Scalar, 6, 6> change_tf_covariance_frame(const typename Eigen::Matrix<Scalar, 6, 6>& covariance_A,
        const typename Eigen::Transform<Scalar, 3, Eigen::Isometry>& transform_B_A) {
    return change_tf_covariance_frame(covariance_A, transform_adjoint(transform_B_A));
}

template<typename Scalar>
inline Eigen::Matrix<Scalar, 6, 1> change_twist_reference_frame(
        const Eigen::Transform<Scalar, 3, Eigen::Isometry>& transform, const Eigen::Matrix<Scalar, 6, 1>& twist) {
    return transform_adjoint(transform) * twist;
}

template<typename Scalar, int D>
inline Eigen::Matrix<Scalar, (D - 1) * 3, (D - 1) * 3> compose_transform_covariance(
        const Eigen::Matrix<Scalar, (D - 1) * 3, (D - 1) * 3>& previous_covariance,
        const Eigen::Matrix<Scalar, (D - 1) * 3, (D - 1) * 3>& relative_covariance,
        const Eigen::Transform<Scalar, D, Eigen::Isometry>& relative_transform,
        const Eigen::Matrix<Scalar, (D - 1) * 3, (D - 1) * 3>& relative_cross_covariance) {
    const Eigen::Matrix<Scalar, (D - 1) * 3, (D - 1) * 3> adj = transform_adjoint(relative_transform.inverse());
    const Eigen::Matrix<Scalar, (D - 1) * 3, (D - 1) * 3> adj_transpose = adj.transpose();
    return adj * previous_covariance * adj_transpose + relative_covariance + relative_cross_covariance * adj_transpose +
           adj * relative_cross_covariance.transpose();
}

template<typename Scalar>
Eigen::Matrix<Scalar, 6, 1> compute_constant_rates(const typename Eigen::Transform<Scalar, 3, Eigen::Isometry>& pose_1,
        const typename Eigen::Transform<Scalar, 3, Eigen::Isometry>& pose_2, const Scalar dt) {
    if (dt <= 0.0) {
        throw std::runtime_error("dt must be > 0.0");
    }
    const typename Eigen::Transform<Scalar, 3, Eigen::Isometry> transform = relative_transform(pose_1, pose_2);
    const typename Eigen::AngleAxis<Scalar> rotation{transform.rotation()};
    return (typename Eigen::Matrix<Scalar, 6, 1>() << rotation.axis() * rotation.angle() / dt,
            transform.translation() / dt)
            .finished();
}

template<typename Scalar, int Dim>
inline Eigen::Transform<Scalar, Dim, Eigen::Isometry> relative_transform(
        const typename Eigen::Transform<Scalar, Dim, Eigen::Isometry>& pose_A_B,
        const typename Eigen::Transform<Scalar, Dim, Eigen::Isometry>& pose_A_C) {
    return pose_A_B.inverse() * pose_A_C;
}

template<typename Derived>
inline Derived rotate_point_covariance(const Eigen::MatrixBase<Derived>& covariance,
        const Eigen::MatrixBase<Derived>& rotation) {
    return rotation * covariance * rotation.transpose();
}

template<typename Derived>
inline Eigen::Matrix<typename Derived::Scalar, 3, 3> rotate_point_covariance(
        const Eigen::Matrix<typename Derived::Scalar, 3, 3>& covariance,
        const Eigen::RotationBase<Derived, 3>& rotation) {
    return rotate_point_covariance(covariance, rotation.toRotationMatrix());
}

template<typename Scalar>
inline Eigen::Matrix<Scalar, 6, 6> transform_adjoint(const Eigen::Transform<Scalar, 3, Eigen::Isometry>& transform) {
    const Eigen::Matrix<Scalar, 3, 1> t = transform.translation();
    const Eigen::Matrix<Scalar, 3, 3> R = transform.rotation();
    return (Eigen::Matrix<Scalar, 6, 6>() << R, Eigen::Matrix<Scalar, 3, 3>::Zero(), skew_symmetric_cross(t) * R, R)
            .finished();
}

}

#endif
