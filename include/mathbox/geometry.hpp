#ifndef MATHBOX_GEOMETRY_HPP
#define MATHBOX_GEOMETRY_HPP

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <mathbox/traits.hpp>
#include <mathbox/types.hpp>

namespace math {

/**
 * @brief Given a relative transform between any two timestamps in the reference frame of A, and a rigid body transform
 * from some other fixed frame B on the body to A, this function computes the relative transform in the reference frame
 * of B.
 *
 * \f[
 *      T_{B1}^{B2} = T_{B1}^{A1} T_{A1}^{A2} T_{A2}^{B2} =  T_B^A T_{A1}^{A2} T_A^B
 * \f]
 *
 * @tparam Scalar
 * @param relative_transform_A \f$T_{A1}^{A2}\f$
 * @param rigid_transform_B_A \f$T_B^A\f$
 * @return Eigen::Transform<Scalar, 3, Eigen::Isometry> \f$T_{B1}^{B2}\f$
 */
template<typename Scalar>
Eigen::Transform<Scalar, 3, Eigen::Isometry> change_relative_transform_frame(
        const typename Eigen::Transform<Scalar, 3, Eigen::Isometry>& relative_transform_A,
        const typename Eigen::Transform<Scalar, 3, Eigen::Isometry>& rigid_transform_B_A);

/**
 * @brief Same as `change_relative_transform_frame`, except \f$T_A^B\f$ is also provided.
 *
 * @tparam Scalar
 * @param relative_transform_A \f$T_{A1}^{A2}\f$
 * @param rigid_transform_B_A \f$T_B^A\f$
 * @param rigid_transform_A_B \f$T_A^B\f$
 * @return Eigen::Transform<Scalar, 3, Eigen::Isometry>
 */
template<typename Scalar>
Eigen::Transform<Scalar, 3, Eigen::Isometry> change_relative_transform_frame(
        const typename Eigen::Transform<Scalar, 3, Eigen::Isometry>& relative_transform_A,
        const typename Eigen::Transform<Scalar, 3, Eigen::Isometry>& rigid_transform_B_A,
        const typename Eigen::Transform<Scalar, 3, Eigen::Isometry>& rigid_transform_A_B);

/**
 * @brief Change the reference frame of a transform covariance matrix from frame A to frame B as
 *
 * \f[
 *      \boldsymbol{\Sigma}_B = Adj_{T_B^A} \boldsymbol{\Sigma}_A (Adj_{T_B^A})^T
 * \f]
 *
 * @tparam Scalar
 * @param covariance_A \f$\boldsymbol{\Sigma}_A\f$
 * @param transform_adjoint_B_A \f$Adj_{T_B^A}\f$
 * @return Eigen::Matrix<Scalar, 6, 6> \f$\boldsymbol{\Sigma}_B\f$
 */
template<typename Scalar>
Eigen::Matrix<Scalar, 6, 6> change_tf_covariance_frame(const typename Eigen::Matrix<Scalar, 6, 6>& covariance_A,
        const typename Eigen::Matrix<Scalar, 6, 6>& transform_adjoint_B_A);

/**
 * @brief Change the reference frame of a transform covariance matrix from frame A to frame B.
 *
 * Note: the function computes the adjoint of the transform, so the other function overload should be preferred if
 * using this function multiple times with the same transform.
 *
 * @tparam Scalar
 * @param covariance_A \f$\boldsymbol{\Sigma}_A\f$
 * @param transform_B_A \f$T_B^A\f$
 * @param translation_before_rotation true if ordering of covariance matrix has translation before rotation
 * @return Eigen::Matrix<Scalar, 6, 6> \f$\boldsymbol{\Sigma}_B\f$
 */
template<typename Scalar>
Eigen::Matrix<Scalar, 6, 6> change_tf_covariance_frame(const typename Eigen::Matrix<Scalar, 6, 6>& covariance_A,
        const typename Eigen::Transform<Scalar, 3, Eigen::Isometry>& transform_B_A,
        const bool translation_before_rotation);

/**
 * @brief Given a transform \f$T_A^B\f$ (i.e. frame B w.r.t. frame A) and the twist in reference frame B, compute
 * the twist in reference frame A. Internally, this function calculates and applies the pose adjoint.
 *
 * \f[
 *      V_A = Ad_{T_A^B} V_B
 * \f]
 *
 * @tparam Scalar
 * @param transform \f$T_A^B\f$
 * @param twist
 * @param translation_before_rotation true if ordering of twist has linear before angular velocity
 * @return Eigen::Matrix<Scalar, 6, 1>
 */
template<typename Scalar>
Eigen::Matrix<Scalar, 6, 1> change_twist_reference_frame(const Eigen::Transform<Scalar, 3, Eigen::Isometry>& transform,
        const Eigen::Matrix<Scalar, 6, 1>& twist, const bool translation_before_rotation);

/**
 * @brief Approximation of the covariance of transform composition.
 *
 * Given the covariance of a (previous) pose \f$T_A^B\f$, \f$\boldsymbol{\Sigma}_{AB}\f$, and a relative transform
 * \f$T_B^C\f$ from that state to a new pose (\f$T_A^C = T_A^B T_B^C\f$) with covariance \f$\boldsymbol{\Sigma}_{BC}\f$,
 * the covariance of the new state \f$\boldsymbol{\Sigma}_{AC}\f$ is computed. The cross correlation between the two
 * relative transforms \f$\boldsymbol{\Sigma}_{AB,BC}\f$ can also be supplied to achieve a better estimate as described
 * in Mangelson et al (2019). Otherwise, the approximation is equal to Barfoot et al's (2013).
 *
 * The ordering of the covariance components is set by `translation_before_rotation`.
 *
 * References:
 *  - Characterizing the Uncertainty of Jointly Distributed Poses in the Lie Algebra, Mangelson et al (2019)
 *  - Associating Uncertainty With Three-Dimensional Poses for Use in Estimation Problems, Barfoot et al (2013)
 *  - https://gtsam.org/2021/02/23/uncertainties-part3.html Eq (15)
 *  - Modern Robotics, Chapter 3, Park & Lynch
 *
 * @tparam Derived
 * @param previous_covariance \f$\boldsymbol{\Sigma}_{AB}\f$
 * @param relative_covariance \f$\boldsymbol{\Sigma}_{BC}\f$
 * @param relative_transform \f$T_B^C\f$
 * @param relative_cross_covariance \f$\boldsymbol{\Sigma}_{AB,BC}\f$
 * @param translation_before_rotation true if ordering of covariance matrices has translation before rotation
 * @return Derived
 */
template<typename Scalar, int D>
    requires(math::is_2d_or_3d<D>)
Eigen::Matrix<Scalar, (D - 1) * 3, (D - 1) * 3> compose_transform_covariance(
        const Eigen::Matrix<Scalar, (D - 1) * 3, (D - 1) * 3>& previous_covariance,
        const Eigen::Matrix<Scalar, (D - 1) * 3, (D - 1) * 3>& relative_covariance,
        const Eigen::Transform<Scalar, D, Eigen::Isometry>& relative_transform,
        const Eigen::Matrix<Scalar, (D - 1) * 3, (D - 1) * 3>& relative_cross_covariance,
        const bool translation_before_rotation);

/**
 * @brief Compute the constant twist (angular and linear velocities) required to pose \f$T_A^B\f$ to $\f$T_A^C\f$ in
 * time `dt`, relative to \f$T_A^B\f$ in frame B (i.e., "body-frame" velocities).
 *
 * Throws error if `dt` is not strictly positive.
 *
 * @tparam Scalar
 * @param pose_1 \f$T_A^B\f$
 * @param pose_2 $\f$T_A^C\f$
 * @param dt
 * @return Eigen::Matrix<Scalar, 6, 1>
 */
template<typename Scalar>
Eigen::Matrix<Scalar, 6, 1> compute_constant_rates(const typename Eigen::Transform<Scalar, 3, Eigen::Isometry>& pose_1,
        const typename Eigen::Transform<Scalar, 3, Eigen::Isometry>& pose_2, const Scalar dt);

/**
 * @brief Compute the general linear interpolation for an SE(3) pose/transform.
 *
 * @tparam Scalar
 * @param T_0
 * @param T_1
 * @param alpha
 * @return Eigen::Transform<Scalar, 3, Eigen::Isometry>
 */
template<typename Scalar>
Eigen::Transform<Scalar, 3, Eigen::Isometry> glerp(const Eigen::Transform<Scalar, 3, Eigen::Isometry>& T_0,
        const Eigen::Transform<Scalar, 3, Eigen::Isometry>& T_1, const Scalar alpha);

/**
 * @brief Compute the relative transform \f$T_B^C = (T_A^B)^{-1} T_A^C\f$ between two poses \f$T_A^B\f$ and \f$T_A^C\f$.
 *
 * @tparam Scalar
 * @tparam D
 * @param pose_A_B
 * @param pose_A_C
 * @return Eigen::Transform<Scalar, D, Eigen::Isometry>
 */
template<typename Scalar, int D>
    requires(math::is_2d_or_3d<D>)
Eigen::Transform<Scalar, D, Eigen::Isometry> relative_transform(
        const typename Eigen::Transform<Scalar, D, Eigen::Isometry>& pose_A_B,
        const typename Eigen::Transform<Scalar, D, Eigen::Isometry>& pose_A_C);

/**
 * @brief Rotate the covariance of an \f$\mathbb{R}^3\f$ point given a rotation matrix.
 *
 * Since points are simply vectors, the application of a rotation matrix means that the returned matrix is
 * \f$R \Sigma R^T\f$, derived from the definition of variance.
 *
 * @tparam Derived
 * @param covariance
 * @param rotation
 * @return Derived
 */
template<typename Derived>
Derived rotate_point_covariance(const Eigen::MatrixBase<Derived>& covariance,
        const Eigen::MatrixBase<Derived>& rotation);

/**
 * @brief Overload of `rotate_point_covariance` for other rotation types.
 *
 * @tparam Derived
 * @param covariance
 * @param rotation
 * @return Eigen::Matrix<typename Derived::Scalar, 3, 3>
 */
template<typename Derived>
Eigen::Matrix<typename Derived::Scalar, 3, 3> rotate_point_covariance(
        const Eigen::Matrix<typename Derived::Scalar, 3, 3>& covariance,
        const Eigen::RotationBase<Derived, 3>& rotation);

template<typename DerivedOmega, typename DerivedV>
    requires(DerivedOmega::RowsAtCompileTime == 1 && DerivedOmega::ColsAtCompileTime == 1 &&
             DerivedV::RowsAtCompileTime == 2 && DerivedV::ColsAtCompileTime == 1 &&
             std::is_same_v<typename DerivedOmega::Scalar, typename DerivedV::Scalar>)
constexpr Eigen::Vector<typename DerivedOmega::Scalar, 2> so_cross(const Eigen::MatrixBase<DerivedOmega>& w,
        const Eigen::MatrixBase<DerivedV>& v);

template<typename DerivedOmega, typename DerivedV>
    requires(DerivedOmega::RowsAtCompileTime == 3 && DerivedOmega::ColsAtCompileTime == 1 &&
             DerivedV::RowsAtCompileTime == 3 && DerivedV::ColsAtCompileTime == 1 &&
             std::is_same_v<typename DerivedOmega::Scalar, typename DerivedV::Scalar>)
constexpr Eigen::Vector<typename DerivedOmega::Scalar, 3> so_cross(const Eigen::MatrixBase<DerivedOmega>& w,
        const Eigen::MatrixBase<DerivedV>& v);

template<typename Derived>
    requires(Derived::RowsAtCompileTime == 2 && Derived::ColsAtCompileTime == 2)
constexpr Eigen::Vector<typename Derived::Scalar, 1> so_from_skew(const Eigen::MatrixBase<Derived>& skew);

template<typename Derived>
    requires(Derived::RowsAtCompileTime == 3 && Derived::ColsAtCompileTime == 3)
constexpr Eigen::Vector<typename Derived::Scalar, 3> so_from_skew(const Eigen::MatrixBase<Derived>& skew);

template<typename Derived>
    requires(Derived::RowsAtCompileTime == 1 && Derived::ColsAtCompileTime == 1)
constexpr Eigen::Matrix<typename Derived::Scalar, 2, 2> so_skew(const Eigen::MatrixBase<Derived>& v);

template<typename Derived>
    requires(Derived::RowsAtCompileTime == 3 && Derived::ColsAtCompileTime == 1)
constexpr Eigen::Matrix<typename Derived::Scalar, 3, 3> so_skew(const Eigen::MatrixBase<Derived>& v);

/**
 * @brief Convert pose from 3D to 2D.
 *
 * @param pose 3D pose to convert
 * @param axis If an axis is zero, it is unused. If it is specified (e.g. Eigen::Vector3d::UnitZ()), then a check is
 * performed that it is equal to the axis of the pose's orientation.
 * @return Eigen::Isometry2d
 */
Pose<2> to_pose_2D(const Pose<3>& pose, const Eigen::Vector3d& axis = Eigen::Vector3d::Zero());

/**
 * @brief Compute the adjoint matrix of a transform T (R, t) in SE(D). It can be used to change the reference frame of
 * twists in the form \f$[v, \omega]\f$ (translation before rotation) or \f$[\omega, v]\f$ (rotation before translation)
 * by \f$V_A = Ad_{T_A^B} V_B\f$.
 *
 * Note that the order matters. If we assume a tangent space where translation is before rotation ([t1, t2, r] for SE(2)
 * or [t1, t2, t3, r1, r2, r3] for SE(3)).
 *
 * References:
 * - Modern Robotics (Lynch & Park)
 *
 * @tparam Scalar
 * @param transform \f$T_A^B\f$
 * @param translation_before_rotation true if ordering is [t, r], false if order is [r, t]
 * @return Eigen::Matrix<Scalar, 6, 6> \f$Ad_{T_A^B}\f$
 */
template<typename Scalar, int D>
    requires(math::is_2d_or_3d<D>)
Eigen::Matrix<Scalar, (D - 1) * 3, (D - 1) * 3> transform_adjoint(
        const Eigen::Transform<Scalar, D, Eigen::Isometry>& transform, const bool translation_before_rotation);

}

#include "mathbox/impl/geometry.hpp"

#endif
