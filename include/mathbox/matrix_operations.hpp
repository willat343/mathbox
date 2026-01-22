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
 * @brief Compute the symmetric matrix \f$ \frac{1}{2} (M + M^T) \f$ of a square matrix \f$ M \f$.
 *
 * This is useful for enforcing symmetry when M should theoretically be symmetric but isn't exactly symmetric due to
 * numerical precision during its construction.
 *
 * @tparam Derived
 * @param m
 * @return constexpr typename Derived::PlainObject
 */
template<typename Derived>
constexpr typename Derived::PlainObject make_symmetric(const Eigen::MatrixBase<Derived>& m);

/**
 * @brief Same as `make_symmetric` but operates in-place on matrix `m`.
 *
 * @tparam Derived
 * @param m
 */
template<typename Derived>
void make_symmetric_inplace(Eigen::MatrixBase<Derived>& m);

/**
 * @brief Return matrix `m` with rows removed where corresponding row of `v` is less than `threshold`. The number of
 * rows in `m` and `v` must be equal.
 *
 * @param m
 * @param v
 * @param threshold
 * @return Eigen::MatrixXd
 */

template<typename DerivedMatrix, typename DerivedVector>
    requires(std::is_same_v<typename DerivedMatrix::Scalar, typename DerivedVector::Scalar>)
Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, DerivedMatrix::ColsAtCompileTime>
remove_rows_by_threshold(const Eigen::MatrixBase<DerivedMatrix>& m, const Eigen::MatrixBase<DerivedVector>& v,
        const typename DerivedMatrix::Scalar threshold);

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
 * @brief Compute the Schur complement matrix H_p and vector b_p for linear system \f$ H_{p} \delta x_k =
 * b_k \f$ from linear system \f$ H \delta x = b\f$:
 * \f[
 *  H_{p} = H_{kk} - H_{km} H_{mm}^{-1} H_{mk}
 * \f]
 * \f[
 *  b_{p} = b_k - H_{km} H_{mm}^{-1} b_m
 * \f]
 *
 * This comes from substituting the first line of the linear system \f$ \delta x_m = H_{mm}^{-1} (b_m - H_{mk} \delta
 * x_k) \f$ into the second line, yielding \f$ (H_{kk} - H_{km} H_{mm}^{-1} H_{mk}) \delta x_k = b_k - H_{km}
 * H_{mm}^{-1} b_m \f$.
 *
 * The input matrix H is assumed symmetric, and output matrix H_p has symmetry enforced.
 *
 * @param H H symmetric matrix
 * @param b b vector
 * @param upper_block_size size of upper block (H_{mm}) or equivalently the index of lower block (H_{kk})
 * @param H_p H_p matrix
 * @param b_p b_p vector
 * @param damping_factor prescaled damping factor to apply to \f$ H_{mm} \f$, typically in [1e-8, 1e-5]
 * @param symmetry_violation_threshold threshold at which numerical symmetry violation is considered an error, compared
 * against the Euclidean norm of the difference of H_p and its transpose
 * @return double the damping applied to H_mm
 */
double schur_complement(const Eigen::Ref<const Eigen::MatrixXd>& H, const Eigen::Ref<const Eigen::VectorXd> b,
        const int upper_block_size, Eigen::MatrixXd& H_p, Eigen::VectorXd& b_p, const double damping_factor,
        const double symmetry_violation_threshold);

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
