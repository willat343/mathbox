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

}

#include "mathbox/impl/matrix_operations.hpp"

#endif
