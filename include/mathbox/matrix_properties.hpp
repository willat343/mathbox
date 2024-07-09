#ifndef MATHBOX_MATRIX_PROPERTIES_HPP
#define MATHBOX_MATRIX_PROPERTIES_HPP

#include <Eigen/Core>
#include <Eigen/Dense>

namespace math {

/**
 * @brief Check if matrix is skew-symmetric.
 *
 * @tparam Derived
 * @param m
 * @return true
 * @return false
 */
template<typename Derived>
bool is_skew_symmetric(const Eigen::DenseBase<Derived>& m);

}

#include "mathbox/impl/matrix_properties.hpp"

#endif
