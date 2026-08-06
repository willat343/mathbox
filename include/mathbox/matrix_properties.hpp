#ifndef MATHBOX_MATRIX_PROPERTIES_HPP
#define MATHBOX_MATRIX_PROPERTIES_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <limits>

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
bool is_positive_definite(const Eigen::MatrixBase<Derived>& m);

/**
 * @brief Check if matrix is positive-semidefinite. This function first checks if the matrix is symmetric, before
 * checking for positive-semidefiniteness through the SelfAdjointEigenSolver.
 *
 * @tparam Derived
 * @param m
 * @return true
 * @return false
 */
template<typename Derived>
bool is_positive_semidefinite(const Eigen::MatrixBase<Derived>& m);

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
        const typename Derived::Scalar precision = std::numeric_limits<typename Derived::Scalar>::epsilon());

/**
 * @brief Compute the relative asymmetry of a square matrix, i.e. \f$ \Vert M - M^T \Vert_F / \Vert M \Vert_F \f$, which
 * is 0 for an exactly symmetric matrix and increases as the matrix deviates from symmetry.
 *
 * Returns 0 if `m` is the zero matrix (which is exactly symmetric), rather than the NaN that would otherwise result
 * from a 0/0 division.
 *
 * @tparam Derived
 * @param m
 * @return typename Derived::Scalar
 */
template<typename Derived>
typename Derived::Scalar relative_asymmetry(const Eigen::MatrixBase<Derived>& m);

/**
 * @brief Check if a matrix is upper triangular, up to some level of precision (default = exact).
 *
 * @tparam Derived
 * @param m
 * @param precision
 * @return true
 * @return false
 */
template<typename Derived>
bool is_upper_triangular(const Eigen::DenseBase<Derived>& m,
        const typename Derived::Scalar precision = std::numeric_limits<typename Derived::Scalar>::epsilon());

/**
 * @brief Count the number of all-zero columns in a matrix, up to some level of precision (default = exact).
 *
 * @tparam Derived
 * @param m
 * @param precision
 * @return int
 */
template<typename Derived>
int num_zero_columns(const Eigen::MatrixBase<Derived>& m,
        const typename Derived::Scalar precision = std::numeric_limits<typename Derived::Scalar>::epsilon());

/**
 * @brief Count the number of all-zero rows in a matrix, up to some level of precision (default = exact).
 *
 * @tparam Derived
 * @param m
 * @param precision
 * @return int
 */
template<typename Derived>
int num_zero_rows(const Eigen::MatrixBase<Derived>& m,
        const typename Derived::Scalar precision = std::numeric_limits<typename Derived::Scalar>::epsilon());

}

#include "mathbox/impl/matrix_properties.hpp"

#endif
