#ifndef MATHBOX_NSPHERE_HPP
#define MATHBOX_NSPHERE_HPP

#include <Eigen/Core>
#include <utility>

namespace math {

/**
 * @brief Convert cartesian coordinates to n-dimensional spherical coordinates, expressed in form \f$\begin{bmatrix}r &
 * \theta_1 & \dots & \theta_{D-1}\end{bmatrix}^T\f$ where \f$r\f$ is the radius, \f$\theta_1, \dots, \theta_{D-2}\f$
 * are angles in range \f$[0, \pi]\f$ and \f$\theta_{D-1}\f$ is in range \f$[-\pi, \pi]\f$, as:
 *
 * \f[
 *      \begin{bmatrix}r \\ \theta_1 \\ \vdots \\ \theta_{D-2} \\ \theta_{D-1}\end{bmatrix} = \begin{bmatrix}
 * \sqrt{x_1^2 + x_2^2 + \dots + x_D^2} \\ \text{atan2}(\sqrt{x_2^2 + \dots + x_D^2}, x_1) \\ \vdots \\
 * \text{atan2}(\sqrt{x_{D-1}^2 + x_D^2}, x_{D-2}) \\ \text{atan2}(x_D, x_{D-1})\end{bmatrix}
 * \f]
 *
 * Note that the `atan2` function here produces angles in the range \f$[-\pi, \pi]\f$.
 *
 * IMPORTANT: This formula is not the same as polar coordinates for D = 2 or spherical coordinates for D = 3. For D = 2,
 * the angular range is usually specified as \f$[0, 2\pi]\f$, not \f$[-\pi, \pi]\f$ as specified here. For D = 3,
 * mathematics often has \f$r, \theta, \phi\f$ representing radius, azimuthal angle (xy-plane) and polar angle (any
 * plane with the z-axis) respectively, and physics often has \f$r, \theta, \phi\f$ representing radius, polar angle
 * (any plane with the z-axis) and azimuthal angle (xy-plane) respectively. One reason for this is how `atan2` is
 * typicall defined. Another reason is because the generalisation to n-dimensions is inconsistent with the mathematics
 * or physics conventions in 3D. If the spherical coordinates are defined on those planes, then the returned cartesian
 * coordinates are \f$\begin{bmatrix}z & x & y\end{bmatrix}^T\f$. Note that the 2D case is consistent with the
 * n-dimensional generalisation, returning \f$\begin{bmatrix}x & y\end{bmatrix}^T\f$.
 *
 * A dimensionality `D` of 1 is a special case, where the spherical coordinates \f$r\f$ and cartesian coordinates
 * \f$x_1\f$ are the same (no square root is performed), and hence the sign disambiguates direction.
 *
 * @tparam D dimensionality (e.g. 1 = line, 2 = circle, 3 = sphere)
 * @tparam Scalar
 * @param cartesian_coordinates
 * @return Eigen::Matrix<Scalar, D, 1>
 */
template<int D, typename Scalar = double>
Eigen::Matrix<Scalar, D, 1> cartesian_to_spherical(const Eigen::Matrix<Scalar, D, 1>& cartesian_coordinates);

/**
 * @brief Same as `cartesian_to_spherical` except only the spherical angles are returned, and range is not.
 *
 * @tparam D dimensionality (e.g. 2 = circle, 3 = sphere)
 * @tparam Scalar
 * @param unit_cartesian_coordinates
 * @return Eigen::Matrix<Scalar, D - 1, 1>
 */
template<int D, typename Scalar = double>
Eigen::Matrix<Scalar, D - 1, 1> cartesian_to_spherical_angles(
        const Eigen::Matrix<Scalar, D, 1>& unit_cartesian_coordinates);

/**
 * @brief Same as `spherical_to_cartesian` except only spherical angles are provided as input, and range is assumed 1.
 *
 * @tparam D dimensionality (e.g. 2 = circle, 3 = sphere)
 * @tparam Scalar
 * @param spherical_angles
 * @return Eigen::Matrix<Scalar, D, 1>
 */
template<int D, typename Scalar = double>
Eigen::Matrix<Scalar, D, 1> spherical_angles_to_unit_cartesian(const Eigen::Matrix<Scalar, D - 1, 1>& spherical_angles);

/**
 * @brief Convert n-dimensional spherical coordinates expressed in form \f$\begin{bmatrix}r & \theta_1 & \dots &
 * \theta_{D-1}\end{bmatrix}^T\f$ where \f$r\f$ is the radius, \f$\theta_1, \dots, \theta_{D-2}\f$ are angles in range
 * \f$[0, \pi]\f$ and \f$theta_{D-1}\f$ is in range\f$[-\pi, \pi]\f$, to cartesian coordinates, as:
 *
 * \f[
 *      \begin{bmatrix}x_1 \\ x_2 \\ \vdots \\ x_{D-1} \\ x_D\end{bmatrix} = r \begin{bmatrix}\cos(\theta_1)
 * \\ \sin(\theta_1) \cos(\theta_2) \\ \vdots \\ \sin(\theta_1) \sin(\theta_2) \dots
 * \sin(\theta_{D-2})\cos(\theta_{D-1}) \\ \sin(\theta_1) \sin(\theta_2) \dots \sin(\theta_{D-2})\sin(\theta_{D-1})
 * \end{bmatrix}
 * \f]
 *
 * IMPORTANT: This formula does not yield the typical spherical coordinates for D = 3 observed in mathematics, which has
 * \f$r, \theta, \phi\f$ representing radius, azimuthal angle (xy-plane) and polar angle (any plane with the z-axis)
 * respectively, nor in physics, which has \f$r, \theta, \phi\f$ representing radius, polar angle (any plane with the
 * z-axis) and azimuthal angle (xy-plane) respectively. This is because the generalisation to n-dimensions is
 * inconsistent with these conventions. If the spherical coordinates are defined on those planes, then the returned
 * cartesian coordinates are \f$\begin{bmatrix}z & x & y\end{bmatrix}^T\f$. Note that the 2D case is consistent with the
 * n-dimensional generalisation, returning \f$\begin{bmatrix}x & y\end{bmatrix}^T\f$.
 *
 * A dimensionality `D` of 1 is a special case, where the spherical coordinates \f$r\f$ and cartesian coordinates
 * \f$x_1\f$ are the same (no square root is performed), and hence the sign disambiguates direction.
 *
 * @tparam D dimensionality (e.g. 1 = line, 2 = circle, 3 = sphere)
 * @tparam Scalar
 * @param spherical_coordinates
 * @return Eigen::Matrix<Scalar, D, 1>
 */
template<int D, typename Scalar = double>
Eigen::Matrix<Scalar, D, 1> spherical_to_cartesian(const Eigen::Matrix<Scalar, D, 1>& spherical_coordinates);

}

#include "mathbox/impl/nsphere.hpp"

#endif
