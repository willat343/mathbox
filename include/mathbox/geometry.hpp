#ifndef MATHBOX_GEOMETRY_HPP
#define MATHBOX_GEOMETRY_HPP

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include <cppbox/enum.hpp>
#include <mathbox/traits.hpp>
#include <mathbox/types.hpp>

namespace math {

/**
 * @brief Degrees of freedom of an orientation in D dimensions, which is the dimension of the tangent space of SO(D).
 *
 * @tparam D
 */
template<int D>
    requires(is_2d_or_3d<D>)
constexpr int SODoF = D * (D - 1) / 2;

/**
 * @brief Degrees of freedom of a pose in D dimensions, which is the dimension of the tangent space of SE(D).
 *
 * @tparam D
 */
template<int D>
    requires(is_2d_or_3d<D>)
constexpr int SEDoF = D + SODoF<D>;

/**
 * @brief Side on which an operation is applied to an element, such as the side of a Lie group element on which a
 * perturbation or another element is composed.
 *
 */
CREATE_SMART_ENUM(ApplicationSide, LEFT, RIGHT)

/**
 * @brief Pose types.
 *
 * JOINT corresponds to SE(n), where the position and orientation are perturbed together.
 * SPLIT corresponds to SO(n) x R^n, where the position and orientation are perturbed separately.
 *
 */
CREATE_SMART_ENUM(SEType, JOINT, SPLIT)

/**
 * @brief Perturbation of a pose \f$ \mathbf{T}_{FC} = (\mathbf{R}, \mathbf{p}) \f$ of a child frame C in a frame F by
 * a tangent vector \f$ \bm{\tau} = (\bm{\rho}, \bm{\theta}) \f$, ordered as position then orientation. The perturbed
 * position and orientation are:
 *  - JOINT, LEFT (\f$ \text{Exp}(\bm{\tau}) \mathbf{T} \f$): \f$ \text{Exp}(\bm{\theta}) \mathbf{p} +
 * \mathbf{V}(\bm{\theta}) \bm{\rho} \f$ and \f$ \text{Exp}(\bm{\theta}) \mathbf{R} \f$
 *  - JOINT, RIGHT (\f$ \mathbf{T} \text{Exp}(\bm{\tau}) \f$): \f$ \mathbf{p} + \mathbf{R} \mathbf{V}(\bm{\theta})
 * \bm{\rho} \f$ and \f$ \mathbf{R} \text{Exp}(\bm{\theta}) \f$
 *  - SPLIT, LEFT: \f$ \mathbf{p} + \bm{\rho} \f$ and \f$ \text{Exp}(\bm{\theta}) \mathbf{R} \f$
 *  - SPLIT, RIGHT: \f$ \mathbf{p} + \bm{\rho} \f$ and \f$ \mathbf{R} \text{Exp}(\bm{\theta}) \f$
 *
 * where \f$ \mathbf{V}(\bm{\theta}) \f$ is the left Jacobian of SO(n). A LEFT orientation perturbation is expressed in
 * F and a RIGHT orientation perturbation is expressed in C. The position perturbation is expressed in C for JOINT,
 * RIGHT and in F otherwise.
 *
 */
struct PosePerturbation {
    SEType type;
    ApplicationSide side;

    bool operator==(const PosePerturbation& other) const;
};

/**
 * @brief Compute the angle between two vectors.
 *
 * @tparam DerivedU
 * @tparam DerivedV
 * @param u
 * @param v
 * @return DerivedU::Scalar angle, in radians, in [0, pi]
 */
template<IsVector DerivedU, IsVector DerivedV>
    requires(std::is_same_v<typename DerivedU::Scalar, typename DerivedV::Scalar> &&
             DerivedU::SizeAtCompileTime == DerivedV::SizeAtCompileTime)
typename DerivedU::Scalar angle_between(const Eigen::MatrixBase<DerivedU>& u, const Eigen::MatrixBase<DerivedV>& v);

/**
 * @brief Convert an angle from degrees to radians.
 *
 * @tparam Scalar
 * @param degrees
 * @return constexpr Scalar
 */
template<typename Scalar>
constexpr Scalar deg2rad(const Scalar degrees);

/**
 * @brief Convert an angle from radians to degrees.
 *
 * @tparam Scalar
 * @param radians
 * @return constexpr Scalar
 */
template<typename Scalar>
constexpr Scalar rad2deg(const Scalar radians);

/**
 * @brief Converts a quaternion to roll-pitch-yaw angles (in radians), in the order of x-axis, y-axis, z-axis (0, 1, 2).
 *
 * @tparam Scalar
 * @param q quaternion
 * @return constexpr Eigen::Vector<Scalar, 3>
 */
template<typename Scalar>
constexpr Eigen::Vector<Scalar, 3> rpy(const Eigen::Quaternion<Scalar>& q);

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
 * @param relative_transform_A \f$ T_{A1}^{A2} \f$
 * @param rigid_transform_B_A \f$ T_B^A \f$
 * @return Eigen::Transform<Scalar, 3, Eigen::Isometry> \f$ T_{B1}^{B2} \f$
 */
template<typename Scalar>
Eigen::Transform<Scalar, 3, Eigen::Isometry> change_relative_transform_frame(
        const typename Eigen::Transform<Scalar, 3, Eigen::Isometry>& relative_transform_A,
        const typename Eigen::Transform<Scalar, 3, Eigen::Isometry>& rigid_transform_B_A);

/**
 * @brief Same as `change_relative_transform_frame`, except \f$ T_A^B \f$ is also provided.
 *
 * @tparam Scalar
 * @param relative_transform_A \f$ T_{A1}^{A2} \f$
 * @param rigid_transform_B_A \f$ T_B^A \f$
 * @param rigid_transform_A_B \f$ T_A^B \f$
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
 *      \boldsymbol{\Sigma}_B = \mathbf{Ad}_{T_B^A} \boldsymbol{\Sigma}_A (\mathbf{Ad}_{T_B^A})^T
 * \f]
 *
 * @tparam Scalar
 * @param covariance_A \f$ \boldsymbol{\Sigma}_A \f$
 * @param adjoint_SE_B_A \f$ \mathbf{Ad}_{T_B^A} \f$
 * @return Eigen::Matrix<Scalar, 6, 6> \f$ \boldsymbol{\Sigma}_B \f$
 */
template<typename Scalar>
Eigen::Matrix<Scalar, 6, 6> change_tf_covariance_frame(const typename Eigen::Matrix<Scalar, 6, 6>& covariance_A,
        const typename Eigen::Matrix<Scalar, 6, 6>& adjoint_SE_B_A);

/**
 * @brief Change the reference frame of a transform covariance matrix from frame A to frame B.
 *
 * This function assumes a translation-rotation ordering in the covariance and adjoint.
 *
 * Note: the function computes the adjoint of the transform, so the other function overload should be preferred if using
 * this function multiple times with the same transform.
 *
 * @tparam Scalar
 * @param covariance_A \f$ \boldsymbol{\Sigma}_A \f$
 * @param transform_B_A \f$ T_B^A \f$
 * @return Eigen::Matrix<Scalar, 6, 6> \f$ \boldsymbol{\Sigma}_B \f$
 */
template<typename Scalar>
Eigen::Matrix<Scalar, 6, 6> change_tf_covariance_frame_tr(const typename Eigen::Matrix<Scalar, 6, 6>& covariance_A,
        const typename Eigen::Transform<Scalar, 3, Eigen::Isometry>& transform_B_A);

/**
 * @brief Overload assuming rotation-translation ordering in the covariance and adjoint.
 *
 * @tparam Scalar
 * @param covariance_A
 * @param transform_B_A
 * @return Eigen::Matrix<Scalar, 6, 6>
 */
template<typename Scalar>
Eigen::Matrix<Scalar, 6, 6> change_tf_covariance_frame_rt(const typename Eigen::Matrix<Scalar, 6, 6>& covariance_A,
        const typename Eigen::Transform<Scalar, 3, Eigen::Isometry>& transform_B_A);

/**
 * @brief Given a transform \f$ T_A^B \f$ (i.e. frame B w.r.t. frame A) and the twist in reference frame B, compute the
 * twist in reference frame A. Internally, this function calculates and applies the pose adjoint.
 *
 * \f[
 *      V_A = Ad_{T_A^B} V_B
 * \f]
 *
 * This function assumes a translation-rotation ordering in the covariance and adjoint.
 *
 * @tparam Scalar
 * @param transform \f$ T_A^B \f$
 * @param twist
 * @param translation_before_rotation true if ordering of twist has linear before angular velocity
 * @return Eigen::Matrix<Scalar, 6, 1>
 */
template<typename Scalar>
Eigen::Matrix<Scalar, 6, 1> change_twist_reference_frame_tr(
        const Eigen::Transform<Scalar, 3, Eigen::Isometry>& transform, const Eigen::Matrix<Scalar, 6, 1>& twist);

/**
 * @brief Overload assuming rotation-translation ordering in the covariance and adjoint.
 *
 * @tparam Scalar
 * @param transform
 * @param twist
 * @return Eigen::Matrix<Scalar, 6, 1>
 */
template<typename Scalar>
Eigen::Matrix<Scalar, 6, 1> change_twist_reference_frame_rt(
        const Eigen::Transform<Scalar, 3, Eigen::Isometry>& transform, const Eigen::Matrix<Scalar, 6, 1>& twist);

/**
 * @brief Approximation of the covariance of transform composition.
 *
 * Given the covariance of a (previous) pose \f$ T_A^B \f$, \f$ \boldsymbol{\Sigma}_{AB} \f$, and a relative transform
 * \f$ T_B^C \f$ from that state to a new pose (\f$ T_A^C = T_A^B T_B^C \f$) with covariance
 * \f$ \boldsymbol{\Sigma}_{BC} \f$, the covariance of the new state \f$ \boldsymbol{\Sigma}_{AC} \f$ is computed. The
 * cross correlation between the two relative transforms \f$ \boldsymbol{\Sigma}_{AB,BC} \f$ can also be supplied to
 * achieve a better estimate as described in Mangelson et al (2019). Otherwise, the approximation is equal to Barfoot et
 * al's (2013).
 *
 * This function assumes a translation-rotation ordering in the covariance and adjoint.
 *
 * References:
 *  - Characterizing the Uncertainty of Jointly Distributed Poses in the Lie Algebra, Mangelson et al (2019)
 *  - Associating Uncertainty With Three-Dimensional Poses for Use in Estimation Problems, Barfoot et al (2013)
 *  - https://gtsam.org/2021/02/23/uncertainties-part3.html Eq (15)
 *  - Modern Robotics, Chapter 3, Park & Lynch
 *
 * @tparam Scalar
 * @tparam D
 * @param previous_covariance \f$ \boldsymbol{\Sigma}_{AB} \f$
 * @param relative_covariance \f$ \boldsymbol{\Sigma}_{BC} \f$
 * @param relative_transform \f$ T_B^C \f$
 * @param relative_cross_covariance \f$ \boldsymbol{\Sigma}_{AB,BC} \f$
 * @return Derived
 */
template<typename Scalar, int D>
    requires(math::is_2d_or_3d<D>)
Eigen::Matrix<Scalar, SEDoF<D>, SEDoF<D>> compose_transform_covariance_tr(
        const Eigen::Matrix<Scalar, SEDoF<D>, SEDoF<D>>& previous_covariance,
        const Eigen::Matrix<Scalar, SEDoF<D>, SEDoF<D>>& relative_covariance,
        const Eigen::Transform<Scalar, D, Eigen::Isometry>& relative_transform,
        const Eigen::Matrix<Scalar, SEDoF<D>, SEDoF<D>>& relative_cross_covariance);

/**
 * @brief Overload assuming rotation-translation ordering in the covariance and adjoint.
 *
 * @tparam Scalar
 * @tparam D
 * @param previous_covariance \f$ \boldsymbol{\Sigma}_{AB} \f$
 * @param relative_covariance \f$ \boldsymbol{\Sigma}_{BC} \f$
 * @param relative_transform \f$ T_B^C \f$
 * @param relative_cross_covariance \f$ \boldsymbol{\Sigma}_{AB,BC} \f$
 * @return Derived
 */
template<typename Scalar, int D>
    requires(math::is_2d_or_3d<D>)
Eigen::Matrix<Scalar, SEDoF<D>, SEDoF<D>> compose_transform_covariance_rt(
        const Eigen::Matrix<Scalar, SEDoF<D>, SEDoF<D>>& previous_covariance,
        const Eigen::Matrix<Scalar, SEDoF<D>, SEDoF<D>>& relative_covariance,
        const Eigen::Transform<Scalar, D, Eigen::Isometry>& relative_transform,
        const Eigen::Matrix<Scalar, SEDoF<D>, SEDoF<D>>& relative_cross_covariance);

/**
 * @brief Compute the constant twist (angular and linear velocities) required to pose \f$ T_A^B \f$ to \f$ T_A^C \f$ in
 * time `dt`, relative to \f$ T_A^B \f$ in frame B (i.e., "body-frame" velocities).
 *
 * Throws error if `dt` is not strictly positive.
 *
 * @tparam Scalar
 * @param pose_1 \f$ T_A^B \f$
 * @param pose_2 \f$ T_A^C \f$
 * @param dt
 * @return Eigen::Matrix<Scalar, 6, 1>
 */
template<typename Scalar>
Eigen::Matrix<Scalar, 6, 1> compute_constant_rates(const typename Eigen::Transform<Scalar, 3, Eigen::Isometry>& pose_1,
        const typename Eigen::Transform<Scalar, 3, Eigen::Isometry>& pose_2, const Scalar dt);

/**
 * @brief Construct a new frame M from two directions given in one frame A, and compute its orientation
 * \f$ \mathbf{R}_A^M \f$ in A. M is the right-handed orthonormal frame whose z axis is the primary direction and whose
 * x axis is the component of the secondary direction orthogonal to it (e.g. a frame with z up and x along a heading):
 * \f[
 *      \mathbf{R}_A^M = [\hat{\mathbf{s}}_\perp, \ \hat{\mathbf{p}} \times \hat{\mathbf{s}}_\perp, \ \hat{\mathbf{p}}],
 *      \qquad \mathbf{s}_\perp = \mathbf{s} - (\mathbf{s} \cdot \hat{\mathbf{p}}) \hat{\mathbf{p}}
 * \f]
 * so the secondary direction is in the xz half-plane of positive x of M.
 *
 * For the rotation between two existing frames in both of which the directions are known, see
 * `rotation_from_directions`. The axes of M are those of the triad \f$ [\hat{\mathbf{p}}, \ \hat{\mathbf{p}} \times
 * \hat{\mathbf{s}}, \ \hat{\mathbf{p}} \times (\hat{\mathbf{p}} \times \hat{\mathbf{s}})] \f$ (normalised) of the TRIAD
 * method referenced there, reordered with the last negated.
 *
 * @tparam Scalar
 * @param primary primary direction \f$ \mathbf{p} \f$, of any non-zero length
 * @param secondary secondary direction \f$ \mathbf{s} \f$, not parallel to the primary direction
 * @return Eigen::Matrix<Scalar, 3, 3> rotation matrix \f$ \mathbf{R}_A^M \f$
 * @throws std::invalid_argument if a direction is zero, or the secondary direction is parallel to the primary
 * direction (to the dummy precision of the scalar)
 */
template<typename Scalar>
Eigen::Matrix<Scalar, 3, 3> frame_from_directions(const Eigen::Vector<Scalar, 3>& primary,
        const Eigen::Vector<Scalar, 3>& secondary);

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
 * @brief Compute the Jacobian \f$ \mathbf{J} = \partial \bm{\tau}_a / \partial \bm{\tau}_b \f$ of the tangent vector
 * \f$ \bm{\tau}_a \f$ of one perturbation of a pose with respect to the tangent vector \f$ \bm{\tau}_b \f$ of another
 * perturbation of the same pose \f$ \mathbf{T}_{FC} = (\mathbf{R}, \mathbf{p}) \f$, where both perturbations result
 * in the same perturbed pose. Tangent vectors \f$ \bm{\tau} = (\bm{\rho}, \bm{\theta}) \f$ are ordered as position
 * then orientation. The Jacobian converts small perturbations and their covariances between conventions:
 *
 * \f[
 *      \bm{\tau}_a \approx \mathbf{J} \bm{\tau}_b, \qquad
 *      \bm{\Sigma}_a \approx \mathbf{J} \bm{\Sigma}_b \mathbf{J}^T
 * \f]
 *
 * It is evaluated at zero perturbation, so these are exact to first order in the perturbation.
 *
 * To first order, where \f$ \mathbf{V}(\bm{\theta}) \approx \mathbf{I} \f$ and \f$ \text{Exp}(\bm{\theta}) \approx
 * \mathbf{I} + [\bm{\theta}]_\times \f$, each perturbation (see `PosePerturbation`) changes the position by
 * \f$ \delta\mathbf{p} \f$ and the orientation by \f$ \mathbf{R} \text{Exp}(\delta\bm{\theta}) \f$, with
 * \f$ (\delta\mathbf{p}, \delta\bm{\theta}) = \mathbf{A} \bm{\tau} \f$ the tangent vector of the equivalent SPLIT,
 * RIGHT perturbation:
 *
 * \f[
 *      \mathbf{A}_\text{SPLIT,RIGHT} = \begin{bmatrix} \mathbf{I} & \mathbf{0} \\ \mathbf{0} & \mathbf{I}
 * \end{bmatrix},
 *      \quad
 *      \mathbf{A}_\text{SPLIT,LEFT} = \begin{bmatrix} \mathbf{I} & \mathbf{0} \\ \mathbf{0} & \mathbf{R}^T
 *      \end{bmatrix}, \quad
 *      \mathbf{A}_\text{JOINT,RIGHT} = \begin{bmatrix} \mathbf{R} & \mathbf{0} \\ \mathbf{0} & \mathbf{I}
 *      \end{bmatrix}, \quad
 *      \mathbf{A}_\text{JOINT,LEFT} = \begin{bmatrix} \mathbf{I} & -[\mathbf{p}]_\times \\ \mathbf{0} & \mathbf{R}^T
 *      \end{bmatrix}
 * \f]
 *
 * These follow from \f$ \text{Exp}(\bm{\theta}) \mathbf{R} = \mathbf{R} \text{Exp}(\mathbf{R}^T \bm{\theta}) \f$ for
 * the LEFT orientation perturbations, \f$ \delta\mathbf{p} = \mathbf{R} \bm{\rho} \f$ for JOINT, RIGHT, and
 * \f$ \delta\mathbf{p} = \bm{\rho} + \bm{\theta} \times \mathbf{p} \f$ for JOINT, LEFT. Equating the equivalent
 * perturbations \f$ \mathbf{A}_a \bm{\tau}_a = \mathbf{A}_b \bm{\tau}_b \f$ gives the Jacobian:
 *
 * \f[
 *      \mathbf{J} = \mathbf{A}_a^{-1} \mathbf{A}_b
 * \f]
 *
 * For example, the Jacobian of JOINT, LEFT with respect to JOINT, RIGHT is the adjoint of the pose (`adjoint_SE_tr`).
 *
 * In 2D the orientation tangent is a scalar \f$ \theta \f$, for which \f$ \mathbf{R}^T \f$ is replaced by 1 and
 * \f$ -[\mathbf{p}]_\times \f$ by \f$ (-p_y, p_x)^T \f$.
 *
 * @tparam Scalar
 * @tparam D
 * @param pose pose \f$ \mathbf{T}_{FC} \f$
 * @param perturbation perturbation \f$ a \f$ of the numerator tangent vector
 * @param with_respect_to perturbation \f$ b \f$ of the denominator tangent vector
 * @return Eigen::Matrix<Scalar, SEDoF<D>, SEDoF<D>> \f$ \mathbf{J} \f$
 */
template<typename Scalar, int D>
    requires(math::is_2d_or_3d<D>)
Eigen::Matrix<Scalar, SEDoF<D>, SEDoF<D>> pose_perturbation_jacobian(
        const Eigen::Transform<Scalar, D, Eigen::Isometry>& pose, const PosePerturbation& perturbation,
        const PosePerturbation& with_respect_to);

/**
 * @brief Compute the rotation vector (the axis scaled by the angle) of the shortest rotation that a unit quaternion
 * represents, \f$ \boldsymbol{\theta} = \mathrm{Log}(\mathbf{q}) \f$, with an angle of at most \f$ \pi \f$. Its
 * inverse, for angles below \f$ \pi \f$, is `rotation_vector_to_quaternion`.
 *
 * With \f$ \mathbf{q} = (w, \mathbf{v}) \f$, sign-flipped if necessary so that \f$ w \geq 0 \f$, this is the
 * angle-axis form (`Eigen::AngleAxis`):
 * \f[
 *      \boldsymbol{\theta} = \frac{2 \operatorname{atan2}(|\mathbf{v}|, w)}{|\mathbf{v}|} \mathbf{v}
 * \f]
 * When \f$ |\mathbf{v}|^2 \f$ is below the machine epsilon, this equals \f$ 2 \mathbf{v} / w \f$ to machine precision,
 * which is used there instead. At zero rotation, the angle-axis form is 0/0, and `Eigen::AngleAxis` returns a zero
 * angle about a fixed axis, whose derivatives are zero, whereas the value and derivatives of \f$ 2 \mathbf{v} / w \f$
 * are exact, e.g. for automatic differentiation.
 *
 * @tparam Scalar scalar type, which may be an automatic differentiation type
 * @param rotation unit quaternion \f$ \mathbf{q} \f$
 * @return Eigen::Vector<Scalar, 3>
 */
template<typename Scalar>
Eigen::Vector<Scalar, 3> quaternion_to_rotation_vector(const Eigen::Quaternion<Scalar>& rotation);

/**
 * @brief Compute the relative transform \f$ T_B^C = (T_A^B)^{-1} T_A^C \f$ between two poses \f$ T_A^B \f$ and
 * \f$ T_A^C \f$.
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
 * @brief Rotate the covariance of an \f$ \mathbb{R}^3 \f$ point given a rotation matrix.
 *
 * Since points are simply vectors, the application of a rotation matrix means that the returned matrix is \f$ R \Sigma
 * R^T \f$, derived from the definition of variance.
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

/**
 * @brief Compute the rotation \f$ \mathbf{R}_A^B \f$ between two existing frames A and B from a primary and a secondary
 * direction that are each known in both frames (the TRIAD method), e.g. the orientation of a body from gravity and a
 * heading known in the body and world frames. The primary direction is matched exactly, and the secondary direction
 * about it:
 * \f[
 *      \mathbf{R}_A^B \hat{\mathbf{p}}_B = \hat{\mathbf{p}}_A, \qquad
 *      \mathbf{R}_A^B \hat{\mathbf{s}}_{B\perp} = \hat{\mathbf{s}}_{A\perp}, \qquad
 *      \mathbf{s}_\perp = \mathbf{s} - (\mathbf{s} \cdot \hat{\mathbf{p}}) \hat{\mathbf{p}}
 * \f]
 * so the angle between the primary and secondary directions may differ between the frames (e.g. from measurement
 * error), and only the component of the secondary direction orthogonal to the primary one is used. The rotation is
 * \f$ \mathbf{R}_A^B = \mathbf{R}_A^M (\mathbf{R}_B^M)^T \f$, with M the frame that the directions define (see
 * `frame_from_directions`), constructed from them in A and in B.
 *
 * References:
 *  - A Passive System for Determining the Attitude of a Satellite, Black (1964), AIAA Journal 2(7), pp. 1350-1351
 *  - Three-Axis Attitude Determination from Vector Observations, Shuster & Oh (1981), Journal of Guidance and Control
 * 4(1), pp. 70-77
 *
 * @tparam Scalar
 * @param primary_A primary direction \f$ \mathbf{p}_A \f$ in frame A, of any non-zero length
 * @param secondary_A secondary direction \f$ \mathbf{s}_A \f$ in frame A, not parallel to the primary direction
 * @param primary_B primary direction \f$ \mathbf{p}_B \f$ in frame B, of any non-zero length
 * @param secondary_B secondary direction \f$ \mathbf{s}_B \f$ in frame B, not parallel to the primary direction
 * @return Eigen::Matrix<Scalar, 3, 3> rotation matrix \f$ \mathbf{R}_A^B \f$
 * @throws std::invalid_argument if a primary direction is zero, or a secondary direction is zero or parallel to its
 * primary direction (to the dummy precision of the scalar)
 */
template<typename Scalar>
Eigen::Matrix<Scalar, 3, 3> rotation_from_directions(const Eigen::Vector<Scalar, 3>& primary_A,
        const Eigen::Vector<Scalar, 3>& secondary_A, const Eigen::Vector<Scalar, 3>& primary_B,
        const Eigen::Vector<Scalar, 3>& secondary_B);

/**
 * @brief Compute the unit quaternion of the rotation that a rotation vector (the axis scaled by the angle)
 * represents, \f$ \mathbf{q} = \mathrm{Exp}(\boldsymbol{\theta}) \f$, the inverse of `quaternion_to_rotation_vector`
 * for angles below \f$ \pi \f$.
 *
 * This is the angle-axis form (`Eigen::AngleAxis`), with the angle \f$ \theta = |\boldsymbol{\theta}| \f$:
 * \f[
 *      \mathbf{q} = \left( \cos\frac{\theta}{2}, \ \sin\frac{\theta}{2} \frac{\boldsymbol{\theta}}{\theta} \right)
 * \f]
 * When \f$ \theta^2 \f$ is below the machine epsilon, this equals \f$ (1, \boldsymbol{\theta} / 2) \f$ to machine
 * precision, which is used there instead. At zero rotation, the angle-axis form is 0/0, whereas the value and
 * derivatives of \f$ (1, \boldsymbol{\theta} / 2) \f$ are exact, e.g. for automatic differentiation.
 *
 * @tparam Scalar scalar type, which may be an automatic differentiation type
 * @param rotation_vector rotation vector \f$ \boldsymbol{\theta} \f$
 * @return Eigen::Quaternion<Scalar>
 */
template<typename Scalar>
Eigen::Quaternion<Scalar> rotation_vector_to_quaternion(const Eigen::Vector<Scalar, 3>& rotation_vector);

/**
 * @brief Compute the SO(2) cross product between an angular velocity `w` and a vector `v`.
 *
 * @tparam DerivedOmega
 * @tparam DerivedV
 * @param w
 * @param v
 * @return Eigen::Vector<typename DerivedOmega::Scalar, 2>
 */
template<typename DerivedOmega, typename DerivedV>
    requires(DerivedOmega::RowsAtCompileTime == 1 && DerivedOmega::ColsAtCompileTime == 1 &&
             DerivedV::RowsAtCompileTime == 2 && DerivedV::ColsAtCompileTime == 1 &&
             std::is_same_v<typename DerivedOmega::Scalar, typename DerivedV::Scalar>)
constexpr Eigen::Vector<typename DerivedOmega::Scalar, 2> so_cross(const Eigen::MatrixBase<DerivedOmega>& w,
        const Eigen::MatrixBase<DerivedV>& v);

/**
 * @brief Compute the SO(3) cross product between an angular velocity `w` and a vector `v`.
 *
 * @tparam DerivedOmega
 * @tparam DerivedV
 * @param w
 * @param v
 * @return Eigen::Vector<typename DerivedOmega::Scalar, 3>
 */
template<typename DerivedOmega, typename DerivedV>
    requires(DerivedOmega::RowsAtCompileTime == 3 && DerivedOmega::ColsAtCompileTime == 1 &&
             DerivedV::RowsAtCompileTime == 3 && DerivedV::ColsAtCompileTime == 1 &&
             std::is_same_v<typename DerivedOmega::Scalar, typename DerivedV::Scalar>)
constexpr Eigen::Vector<typename DerivedOmega::Scalar, 3> so_cross(const Eigen::MatrixBase<DerivedOmega>& w,
        const Eigen::MatrixBase<DerivedV>& v);

/**
 * @brief Compute the so(2) vector from its skew-symmetric matrix representation.
 *
 * @tparam Derived
 * @param skew
 * @return Eigen::Vector<typename Derived::Scalar, 1>
 */
template<typename Derived>
    requires(Derived::RowsAtCompileTime == 2 && Derived::ColsAtCompileTime == 2)
constexpr Eigen::Vector<typename Derived::Scalar, 1> so_from_skew(const Eigen::MatrixBase<Derived>& skew);

/**
 * @brief Compute the so(3) vector from its skew-symmetric matrix representation.
 *
 * @tparam Derived
 * @param skew
 * @return Eigen::Vector<typename Derived::Scalar, 3>
 */
template<typename Derived>
    requires(Derived::RowsAtCompileTime == 3 && Derived::ColsAtCompileTime == 3)
constexpr Eigen::Vector<typename Derived::Scalar, 3> so_from_skew(const Eigen::MatrixBase<Derived>& skew);

/**
 * @brief Compute the skew symmetric matrix for the 2D so Lie algebra.
 *
 * @tparam Derived
 * @param v
 * @return Eigen::Matrix<typename Derived::Scalar, 2, 2>
 */
template<typename Derived>
    requires(Derived::RowsAtCompileTime == 1 && Derived::ColsAtCompileTime == 1)
constexpr Eigen::Matrix<typename Derived::Scalar, 2, 2> so_skew(const Eigen::MatrixBase<Derived>& v);

/**
 * @brief Compute the skew symmetric matrix for the 3D so Lie algebra.
 *
 * Equivalent to `skew_symmetric_cross(v)`.
 *
 * @tparam Derived
 * @param v
 * @return Eigen::Matrix<typename Derived::Scalar, 3, 3>
 */
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

template<int D>
    requires(is_2d_or_3d<D>)
Pose<D> to_pose_ND(const Pose<3>& pose, [[maybe_unused]] const Eigen::Vector3d& axis = Eigen::Vector3d::Zero());

/**
 * @brief Compute the adjoint matrix of a transform T (R, t) in SE(D). It can be used to change the reference frame of
 * twists in the form \f$ [v, \omega] \f$ (translation before rotation) or \f$ [\omega, v] \f$ (rotation before
 * translation) by \f$ V_A = Ad_{T_A^B} V_B \f$.
 *
 * Note that the order matters. If we assume a tangent space where translation is before rotation ([t1, t2, r] for SE(2)
 * or [t1, t2, t3, r1, r2, r3] for SE(3)).
 *
 * This function assumes a translation-rotation [t, r] ordering.
 *
 * References: - Modern Robotics (Lynch & Park)
 *
 * @tparam Scalar
 * @tparam D
 * @param transform \f$ T_A^B \f$
 * @return Eigen::Matrix<Scalar, 6, 6> \f$ Ad_{T_A^B} \f$
 */
template<typename Scalar, int D>
    requires(math::is_2d_or_3d<D>)
Eigen::Matrix<Scalar, SEDoF<D>, SEDoF<D>> adjoint_SE_tr(const Eigen::Transform<Scalar, D, Eigen::Isometry>& transform);

/**
 * @brief Overload assuming rotation-translation [r, t] ordering.
 *
 * @tparam Scalar
 * @tparam D
 * @param transform \f$ T_A^B \f$
 * @return Eigen::Matrix<Scalar, 6, 6> \f$ Ad_{T_A^B} \f$
 */
template<typename Scalar, int D>
    requires(math::is_2d_or_3d<D>)
Eigen::Matrix<Scalar, SEDoF<D>, SEDoF<D>> adjoint_SE_rt(const Eigen::Transform<Scalar, D, Eigen::Isometry>& transform);

}

#include "mathbox/impl/geometry.hpp"

#endif
