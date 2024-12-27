#ifndef MATHBOX_MATRIX_OPERATIONS_HPP
#define MATHBOX_MATRIX_OPERATIONS_HPP

#include <Eigen/Core>
#include <Eigen/Dense>

namespace math {

/**
 * @brief Compute the cumulative sum of a matrix along its rows from right to left. The largest entry is in the left.
 *
 * Example:
 * \f[
 *      \begin{bmatrix}1 & 2\\ 3 & 4\end{bmatrix} \rightarrow \begin{bmatrix}3 & 2\\ 7 & 4\end{bmatrix}
 * \f]
 *
 * @tparam Derived
 * @param m
 * @return constexpr Derived
 */
template<typename Derived>
constexpr Derived cumulative_row_left_sum(const Eigen::DenseBase<Derived>& m);

/**
 * @brief Compute the cumulative sum of a matrix along its columns from bottom to top. The largest entry is in the top.
 *
 * Example:
 * \f[
 *      \begin{bmatrix}1 & 2\\ 3 & 4\end{bmatrix} \rightarrow \begin{bmatrix}4 & 6\\ 3 & 4\end{bmatrix}
 * \f]
 *
 * @tparam Derived
 * @param m
 * @return constexpr Derived
 */
template<typename Derived>
constexpr Derived cumulative_col_top_sum(const Eigen::DenseBase<Derived>& m);

/**
 * @brief Re-order a symmetric matrix (e.g. covariance matrix) by swapping the blocks according to some boundary index.
 *
 * \f[
 *      \begin{bmatrix}A & X \\ X^T & B\end{bmatrix} \rightarrow \begin{bmatrix}B & X^T \\ X & A\end{bmatrix}
 * \f]
 *
 * Throws error if boundary is 0 or >= matrix size, or if matrix is not square. The symmetry of the matrix is not
 * checked.
 *
 * @tparam Derived
 * @param m matrix
 * @param boundary B block row index
 * @return Derived
 */
template<typename Derived>
Derived reorder_symmetric_matrix(const Eigen::MatrixBase<Derived>& m, const Eigen::Index boundary);

/**
 * @brief Compute the 3D skew-symmetric matrix \f$[\mathbf{x}]_\times\f$ which when multiplied with a vector, is
 * equivalent to the vector cross product.
 *
 * \f[
 *      \mathbf{x} \times \mathbf{y} = [\mathbf{x}]_\times \mathbf{y}
 * \f]
 *
 * The matrix has form
 *
 * \f[
 *      [\mathbf{x}]_\times = \begin{bmatrix}0 & -x_3 & x_2 \\ x_3 & 0 & -x_1 \\ -x_2 & x_1 & 0\end{bmatrix}
 * \f]
 *
 * Note also that \f$\mathbf{x} \times \mathbf{y} = - \mathbf{y} \times \mathbf{x} = - [\mathbf{y}]_\times \mathbf{x}
 * = [\mathbf{y}]_\times^T \mathbf{x} = [\mathbf{x}]_\times \mathbf{y}\f$.
 *
 * @tparam Scalar
 * @param v
 * @return Eigen::Matrix<Scalar, 3, 3>
 */
template<typename Scalar>
Eigen::Matrix<Scalar, 3, 3> skew_symmetric_cross(const Eigen::Matrix<Scalar, 3, 1>& v);

}

#include "mathbox/impl/matrix_operations.hpp"

#endif
