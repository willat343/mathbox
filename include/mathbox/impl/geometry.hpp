#ifndef MATHBOX_IMPL_GEOMETRY_HPP
#define MATHBOX_IMPL_GEOMETRY_HPP

#include <cppbox/exceptions.hpp>

#include "mathbox/geometry.hpp"
#include "mathbox/lerp.hpp"
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
        const typename Eigen::Transform<Scalar, 3, Eigen::Isometry>& transform_B_A,
        const bool translation_before_rotation) {
    return change_tf_covariance_frame(covariance_A, transform_adjoint(transform_B_A, translation_before_rotation));
}

template<typename Scalar>
inline Eigen::Matrix<Scalar, 6, 1> change_twist_reference_frame(
        const Eigen::Transform<Scalar, 3, Eigen::Isometry>& transform, const Eigen::Matrix<Scalar, 6, 1>& twist,
        const bool translation_before_rotation) {
    return transform_adjoint(transform, translation_before_rotation) * twist;
}

template<typename Scalar, int D>
    requires(math::is_2d_or_3d<D>)
inline Eigen::Matrix<Scalar, (D - 1) * 3, (D - 1) * 3> compose_transform_covariance(
        const Eigen::Matrix<Scalar, (D - 1) * 3, (D - 1) * 3>& previous_covariance,
        const Eigen::Matrix<Scalar, (D - 1) * 3, (D - 1) * 3>& relative_covariance,
        const Eigen::Transform<Scalar, D, Eigen::Isometry>& relative_transform,
        const Eigen::Matrix<Scalar, (D - 1) * 3, (D - 1) * 3>& relative_cross_covariance,
        const bool translation_before_rotation) {
    const Eigen::Matrix<Scalar, (D - 1) * 3, (D - 1) * 3> adj =
            transform_adjoint(relative_transform.inverse(), translation_before_rotation);
    const Eigen::Matrix<Scalar, (D - 1) * 3, (D - 1) * 3> adj_transpose = adj.transpose();
    return adj * previous_covariance * adj_transpose + relative_covariance + relative_cross_covariance * adj_transpose +
           adj * relative_cross_covariance.transpose();
}

template<typename Scalar>
Eigen::Matrix<Scalar, 6, 1> compute_constant_rates(const typename Eigen::Transform<Scalar, 3, Eigen::Isometry>& pose_1,
        const typename Eigen::Transform<Scalar, 3, Eigen::Isometry>& pose_2, const Scalar dt) {
    throw_if(dt <= 0.0, "dt must be > 0.0");
    const typename Eigen::Transform<Scalar, 3, Eigen::Isometry> transform = relative_transform(pose_1, pose_2);
    const typename Eigen::AngleAxis<Scalar> rotation{transform.rotation()};
    return (typename Eigen::Matrix<Scalar, 6, 1>() << rotation.axis() * rotation.angle() / dt,
            transform.translation() / dt)
            .finished();
}

template<typename Scalar>
inline Eigen::Transform<Scalar, 3, Eigen::Isometry> glerp(const Eigen::Transform<Scalar, 3, Eigen::Isometry>& T_0,
        const Eigen::Transform<Scalar, 3, Eigen::Isometry>& T_1, const Scalar alpha) {
    const Eigen::Translation<Scalar, 3> t_lerp = Eigen::Translation<Scalar, 3>{lerp(
            Eigen::Matrix<Scalar, 3, 1>{T_0.translation()}, Eigen::Matrix<Scalar, 3, 1>{T_1.translation()}, alpha)};
    const Eigen::Quaternion<Scalar> q_lerp =
            Eigen::Quaternion<Scalar>{T_0.rotation()}.slerp(alpha, Eigen::Quaternion<Scalar>{T_1.rotation()});
    return t_lerp * q_lerp;
}

template<typename Scalar, int D>
    requires(math::is_2d_or_3d<D>)
inline Eigen::Transform<Scalar, D, Eigen::Isometry> relative_transform(
        const typename Eigen::Transform<Scalar, D, Eigen::Isometry>& pose_A_B,
        const typename Eigen::Transform<Scalar, D, Eigen::Isometry>& pose_A_C) {
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

template<typename DerivedOmega, typename DerivedV>
    requires(DerivedOmega::RowsAtCompileTime == 1 && DerivedOmega::ColsAtCompileTime == 1 &&
             DerivedV::RowsAtCompileTime == 2 && DerivedV::ColsAtCompileTime == 1 &&
             std::is_same_v<typename DerivedOmega::Scalar, typename DerivedV::Scalar>)
constexpr Eigen::Vector<typename DerivedOmega::Scalar, 2> so_cross(const Eigen::MatrixBase<DerivedOmega>& w,
        const Eigen::MatrixBase<DerivedV>& v) {
    return w[0] * Eigen::Vector<typename DerivedOmega::Scalar, 2>{-v[1], v[0]};
}

template<typename DerivedOmega, typename DerivedV>
    requires(DerivedOmega::RowsAtCompileTime == 3 && DerivedOmega::ColsAtCompileTime == 1 &&
             DerivedV::RowsAtCompileTime == 3 && DerivedV::ColsAtCompileTime == 1 &&
             std::is_same_v<typename DerivedOmega::Scalar, typename DerivedV::Scalar>)
constexpr Eigen::Vector<typename DerivedOmega::Scalar, 3> so_cross(const Eigen::MatrixBase<DerivedOmega>& w,
        const Eigen::MatrixBase<DerivedV>& v) {
    return w.cross(v);
}

template<typename Derived>
    requires(Derived::RowsAtCompileTime == 2 && Derived::ColsAtCompileTime == 2)
constexpr inline Eigen::Vector<typename Derived::Scalar, 1> so_from_skew(const Eigen::MatrixBase<Derived>& skew) {
    return Eigen::Vector<typename Derived::Scalar, 1>{m(1, 1)};
}

template<typename Derived>
    requires(Derived::RowsAtCompileTime == 3 && Derived::ColsAtCompileTime == 3)
constexpr inline Eigen::Vector<typename Derived::Scalar, 3> so_from_skew(const Eigen::MatrixBase<Derived>& skew) {
    return Eigen::Vector<typename Derived::Scalar, 3>{m(2, 1), m(0, 2), m(1, 0)};
}

template<typename Derived>
    requires(Derived::RowsAtCompileTime == 1 && Derived::ColsAtCompileTime == 1)
constexpr inline Eigen::Matrix<typename Derived::Scalar, 2, 2> so_skew(const Eigen::MatrixBase<Derived>& v) {
    return (Eigen::Matrix<typename Derived::Scalar, 2, 2>() << static_cast<typename Derived::Scalar>(0), -v.value(),
            static_cast<typename Derived::Scalar>(0), v.value())
            .finished();
}

template<typename Derived>
    requires(Derived::RowsAtCompileTime == 3 && Derived::ColsAtCompileTime == 1)
constexpr inline Eigen::Matrix<typename Derived::Scalar, 3, 3> so_skew(const Eigen::MatrixBase<Derived>& v) {
    return (Eigen::Matrix<typename Derived::Scalar, 3, 3>() << static_cast<typename Derived::Scalar>(0), -v[2], v[1],
            v[2], static_cast<typename Derived::Scalar>(0), -v[0], -v[1], v[0],
            static_cast<typename Derived::Scalar>(0))
            .finished();
}

inline Pose<2> to_pose_2D(const Pose<3>& pose, const Eigen::Vector3d& axis) {
    const Eigen::AngleAxisd orientation{pose.rotation()};
    throw_if(!axis.isZero() && orientation.axis().isApprox(axis),
            "Orientation must be purely a rotation about a normalized axis argument (e.g. UnitZ) to be convertible to "
            "Eigen::Rotation2D, or this check can be disabled by passing the default value for axis (zero).");
    return Eigen::Translation2d{pose.translation()[0], pose.translation()[1]} * Eigen::Rotation2Dd{orientation.angle()};
}

template<typename Scalar, int D>
    requires(math::is_2d_or_3d<D>)
inline Eigen::Matrix<Scalar, (D - 1) * 3, (D - 1) * 3> transform_adjoint(
        const Eigen::Transform<Scalar, D, Eigen::Isometry>& transform, const bool translation_before_rotation) {
    const Eigen::Matrix<Scalar, D, 1> t = transform.translation();
    const Eigen::Matrix<Scalar, D, D> R = transform.rotation();
    if constexpr (D == 2) {
        if (translation_before_rotation) {
            return (Eigen::Matrix<Scalar, (D - 1) * 3, (D - 1) * 3>() << R, Eigen::Matrix<Scalar, D, 1>(t[1], -t[0]),
                    Eigen::Matrix<Scalar, 1, D>::Zero(), 1)
                    .finished();
        } else {
            return (Eigen::Matrix<Scalar, (D - 1) * 3, (D - 1) * 3>() << R, Eigen::Matrix<Scalar, D, 1>::Zero(),
                    Eigen::Matrix<Scalar, 1, D>(t[1], -t[0]), 1)
                    .finished();
        }
    } else {
        if (translation_before_rotation) {
            return (Eigen::Matrix<Scalar, (D - 1) * 3, (D - 1) * 3>() << R, skew_symmetric_cross(t) * R,
                    Eigen::Matrix<Scalar, D, D>::Zero(), R)
                    .finished();
        } else {
            return (Eigen::Matrix<Scalar, (D - 1) * 3, (D - 1) * 3>() << R, Eigen::Matrix<Scalar, D, D>::Zero(),
                    skew_symmetric_cross(t) * R, R)
                    .finished();
        }
    }
}

}

#if !MATHBOX_HEADER_ONLY
namespace math {

extern template Eigen::Isometry3d change_relative_transform_frame<double>(const Eigen::Isometry3d&,
        const Eigen::Isometry3d&);

extern template Eigen::Isometry3d change_relative_transform_frame<double>(const Eigen::Isometry3d&,
        const Eigen::Isometry3d&, const Eigen::Isometry3d&);

extern template Eigen::Matrix<double, 6, 6> change_tf_covariance_frame<double>(const Eigen::Matrix<double, 6, 6>&,
        const Eigen::Isometry3d&, const bool);

extern template Eigen::Matrix<double, 6, 1> change_twist_reference_frame<double>(const Eigen::Isometry3d&,
        const Eigen::Matrix<double, 6, 1>&, const bool);

extern template Eigen::Matrix<double, 3, 3> compose_transform_covariance<double, 2>(
        const Eigen::Matrix<double, 3, 3>& previous_covariance, const Eigen::Matrix<double, 3, 3>& relative_covariance,
        const Eigen::Isometry2d& relative_transform, const Eigen::Matrix<double, 3, 3>& relative_cross_covariance,
        const bool);
extern template Eigen::Matrix<double, 6, 6> compose_transform_covariance<double, 3>(
        const Eigen::Matrix<double, 6, 6>& previous_covariance, const Eigen::Matrix<double, 6, 6>& relative_covariance,
        const Eigen::Isometry3d& relative_transform, const Eigen::Matrix<double, 6, 6>& relative_cross_covariance,
        const bool);

extern template Eigen::Matrix<double, 6, 1> compute_constant_rates<double>(const typename Eigen::Isometry3d& pose_1,
        const typename Eigen::Isometry3d& pose_2, const double dt);

extern template Eigen::Isometry3d glerp<double>(const Eigen::Isometry3d& T_0, const Eigen::Isometry3d& T_1,
        const double alpha);

extern template Eigen::Isometry2d relative_transform<double, 2>(const typename Eigen::Isometry2d& pose_A_B,
        const typename Eigen::Isometry2d& pose_A_C);
extern template Eigen::Isometry3d relative_transform<double, 3>(const typename Eigen::Isometry3d& pose_A_B,
        const typename Eigen::Isometry3d& pose_A_C);

extern template Eigen::Matrix2d rotate_point_covariance<Eigen::Matrix2d>(
        const Eigen::MatrixBase<Eigen::Matrix2d>& covariance, const Eigen::MatrixBase<Eigen::Matrix2d>& rotation);
extern template Eigen::Matrix3d rotate_point_covariance<Eigen::Matrix3d>(
        const Eigen::MatrixBase<Eigen::Matrix3d>& covariance, const Eigen::MatrixBase<Eigen::Matrix3d>& rotation);

extern template Eigen::Vector2d so_cross<Eigen::Vector<double, 1>, Eigen::Vector2d>(
        const Eigen::MatrixBase<Eigen::Vector<double, 1>>&, const Eigen::MatrixBase<Eigen::Vector2d>&);
extern template Eigen::Vector3d so_cross<Eigen::Vector3d, Eigen::Vector3d>(const Eigen::MatrixBase<Eigen::Vector3d>&,
        const Eigen::MatrixBase<Eigen::Vector3d>&);

extern template Eigen::Vector<double, 1> so_from_skew<Eigen::Matrix2d>(const Eigen::MatrixBase<Eigen::Matrix2d>&);
extern template Eigen::Vector3d so_from_skew<Eigen::Matrix3d>(const Eigen::MatrixBase<Eigen::Matrix3d>&);

extern template Eigen::Matrix2d so_skew<Eigen::Vector<double, 1>>(const Eigen::MatrixBase<Eigen::Vector<double, 1>>&);
extern template Eigen::Matrix3d so_skew<Eigen::Vector3d>(const Eigen::MatrixBase<Eigen::Vector3d>&);

extern template Eigen::Matrix<double, 3, 3> transform_adjoint<double, 2>(const Eigen::Isometry2d& transform,
        const bool);
extern template Eigen::Matrix<double, 6, 6> transform_adjoint<double, 3>(const Eigen::Isometry3d& transform,
        const bool);

}
#endif

#endif
