#ifndef MATHBOX_VECTOR_OPERATIONS_HPP
#define MATHBOX_VECTOR_OPERATIONS_HPP

#include <Eigen/Core>
#include <Eigen/Dense>

namespace math {

/**
 * @brief Generate a linearly spaced dynamically-sized matrix from an interpolation step. The number of points are
 * chosen such that the actual interpolation step is greater than or equal to `step`.
 *
 * The first element will be `start` and the last element will be `end`.
 *
 * @tparam Scalar
 * @param step
 * @param start
 * @param end
 * @return Eigen::Matrix<Scalar, Eigen::Dynamic, 1>
 */
template<typename Scalar>
Eigen::Matrix<Scalar, Eigen::Dynamic, 1> lin_spaced(const Scalar step, const Scalar start, const Scalar end);

/**
 * @brief Generate a range from `start` to `end` incrementing by `step`.
 *
 * The first element will be `start` and the last element will be `end` minus the modulus of (`end` - `start`) and
 * `step`. It is thus at most `end`.
 *
 * @tparam Scalar
 * @param step
 * @param start
 * @param end
 * @return Eigen::Matrix<Scalar, Eigen::Dynamic, 1>
 */
template<typename Scalar>
Eigen::Matrix<Scalar, Eigen::Dynamic, 1> range(const Scalar step, const Scalar start, const Scalar end);

}

#include "mathbox/impl/vector_operations.hpp"

#endif
