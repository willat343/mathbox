#ifndef MATHBOX_MATRIX_PROPERTIES_HPP
#define MATHBOX_MATRIX_PROPERTIES_HPP

#include <Eigen/Core>
#include <Eigen/Dense>

namespace math {

/**
 * @brief Check if matrix has all diagonals positive (strictly greater than zero).
 * 
 * @tparam Derived 
 * @param m 
 * @return true 
 * @return false 
 */
template<typename Derived>
bool has_positive_diagonals(const Eigen::MatrixBase<Derived>& m);

/**
 * @brief Check if matrix is positive-definite. This function first checks if the matrix is symmetric, before checking
 * for positive-definiteness through Cholesky decomposition.
 *
 * @tparam Derived
 * @param m
 * @return true
 * @return false
 */
template<typename Derived>
bool is_positive_definite(const Eigen::EigenBase<Derived>& m);

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

/**
 * @brief Check if a matrix is symmetric about its diagonal, up to some level of precision (default = exact).
 * 
 * @tparam Derived 
 * @param m 
 * @param precision 
 * @return true 
 * @return false 
 */
template<typename Derived>
bool is_symmetric(const Eigen::DenseBase<Derived>& m,
        const typename Derived::Scalar precision = static_cast<typename Derived::Scalar>(0));

}

#include "mathbox/impl/matrix_properties.hpp"

#endif
