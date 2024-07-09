#ifndef MATHBOX_GEOMETRY_HPP
#define MATHBOX_GEOMETRY_HPP

#include <Eigen/Core>

namespace math {

/**
 * @brief Convert a vector representing an angle-axis rotation in form \f$\mathbf{r} = \Vert \mathbf{r} \Vert
 * \mathbf{\hat{r}}\f$ to angle-axis object. If the angle \f$\Vert \mathbf{r} \Vert == 0\f$, then the unit X axis is
 * used.
 *
 * @tparam Scalar
 * @param vector
 * @return Eigen::AngleAxis<Scalar>
 */
template<typename Scalar>
Eigen::AngleAxis<Scalar> to_angleaxis(const Eigen::Matrix<Scalar, 3, 1>& vector);

/**
 * @brief Convert an angle-axis object to a vector \f$\mathbf{r}\f$ by multiplying the angle \f$\Vert \mathbf{r}
 * \Vert\f$ by the rotation axis \f$\mathbf{\hat{r}}\f$.
 *
 * @tparam Scalar
 * @param angleaxis
 * @return Eigen::Matrix<Scalar, 3, 1>
 */
template<typename Scalar>
Eigen::Matrix<Scalar, 3, 1> to_vector(const Eigen::AngleAxis<Scalar>& angleaxis);

}

#include "mathbox/impl/geometry.hpp"

#endif
