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
 * @brief Return matrix `m` with rows removed where corresponding row of `v` is less than `threshold`. The number of
 * rows in `m` and `v` must be equal.
 *
 * @tparam DerivedMatrix
 * @tparam DerivedVector
 * @param m
 * @param v
 * @param threshold
 * @return Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, DerivedMatrix::ColsAtCompileTime>
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
 * @brief Compute the Schur complement matrix H_p and vector b_p for linear system \f$ H_{p} \delta x_k = b_k \f$ from
 * linear system \f$ H \delta x = b \f$:
 * \f[
 *      H_{p} = H_{kk} - H_{km} H_{mm}^{-1} H_{mk}
 * \f]
 * \f[
 *      b_{p} = b_k - H_{km} H_{mm}^{-1} b_m
 * \f]
 *
 * This comes from substituting the first line of the linear system \f$ \delta x_m = H_{mm}^{-1} (b_m - H_{mk} \delta
 * x_k) \f$ into the second line, yielding \f$ (H_{kk} - H_{km} H_{mm}^{-1} H_{mk}) \delta x_k = b_k - H_{km}
 * H_{mm}^{-1} b_m \f$.
 *
 * The input matrix H is assumed symmetric, and output matrix H_p has symmetry enforced.
 *
 * Damping \f$ \lambda = \text{damping_factor} \cdot \text{tr}(H_{mm}) / \text{rows}(H_{mm}) \f$ is added to the
 * diagonal of \f$ H_{mm} \f$, which its Cholesky (LLT) decomposition requires to be positive definite. It is only
 * needed when \f$ H_{mm} \f$ is singular or numerically indefinite. It is not free:
 * \f$ (H_{mm} + \lambda I)^{-1} \prec H_{mm}^{-1} \f$ under-estimates the subtracted term, so \f$ H_p \f$ is larger
 * than the true Schur complement. Being trace-scaled, it is never small for a large-magnitude \f$ H_{mm} \f$.
 *
 * Jacobi scaling (diagonal equilibration) factorises \f$ \tilde{H}_{mm} = D^{-1} H_{mm} D^{-1} \f$ with
 * \f$ D = \text{diag}(H_{mm})^{1/2} \f$, which has a unit diagonal. Substituting
 * \f$ H_{mm}^{-1} = D^{-1} \tilde{H}_{mm}^{-1} D^{-1} \f$ gives:
 * \f[
 *      H_{p} = H_{kk} - (H_{km} D^{-1}) \tilde{H}_{mm}^{-1} (D^{-1} H_{mk})
 * \f]
 * \f[
 *      b_{p} = b_k - (H_{km} D^{-1}) \tilde{H}_{mm}^{-1} (D^{-1} b_m)
 * \f]
 * Unlike damping this is exact: the \f$ D \f$ factors cancel, \f$ H_{kk} \f$ and \f$ b_k \f$ are untouched, and
 * \f$ H_p \f$ and \f$ b_p \f$ are produced in the original units. It is therefore purely a change of the
 * floating-point path, and is also a no-op for any damping, since scaling the damped block is equivalent to adding
 * \f$ \lambda / H_{mm}(i, i) \f$ to the unit diagonal. It is useful when the variables of \f$ H_{mm} \f$ have widely
 * differing units, as van der Sluis (1969) showed that scaling to a unit diagonal is within a factor of the matrix
 * dimension of the optimal diagonal scaling, \f$ \kappa(DAD) \leq n \min_{D'} \kappa(D'AD') \f$. For a positive
 * semi-definite \f$ H \f$, Cauchy-Schwarz gives \f$ |H_{ij}| \leq \sqrt{H_{ii} H_{jj}} \f$, so the scaled entries
 * satisfy \f$ |\tilde{H}_{ij}| \leq 1 \f$ and cannot overflow.
 *
 * The computation is in the scalar type of H, e.g. long double for a reference of higher precision.
 *
 * @tparam DerivedMatrix
 * @tparam DerivedVector
 * @param H H symmetric matrix
 * @param b b vector
 * @param upper_block_size size of upper block (H_{mm}) or equivalently the index of lower block (H_{kk})
 * @param H_p H_p matrix
 * @param b_p b_p vector
 * @param damping_factor prescaled damping factor to apply to \f$ H_{mm} \f$, recommended to be 0.0 unless the
 * decomposition fails, because it biases \f$ H_p \f$ as described above
 * @param symmetry_violation_threshold threshold at which numerical symmetry violation is considered an error, compared
 * against the relative asymmetry \f$ \| H_p - H_p^T \| / \| H_p \| \f$ (Frobenius norms). It is recommended to be a
 * small non-zero value such as 1.0e-9, since an exactly symmetric H_p is not achievable in floating-point arithmetic,
 * while a violation far above the round-off indicates an asymmetric H
 * @param jacobi_scaling scale \f$ H_{mm} \f$ to a unit diagonal before its decomposition, which is exact and only
 * changes the numerical conditioning of the decomposition. It is recommended to be false, because the Cholesky
 * backward error is componentwise and therefore already invariant to a diagonal scaling, and empirically it has not
 * been measured to improve accuracy.
 * @return typename DerivedMatrix::Scalar the damping applied to H_mm
 */
template<typename DerivedMatrix, typename DerivedVector>
    requires(std::is_same_v<typename DerivedMatrix::Scalar, typename DerivedVector::Scalar>)
typename DerivedMatrix::Scalar schur_complement(const Eigen::MatrixBase<DerivedMatrix>& H,
        const Eigen::MatrixBase<DerivedVector>& b, const int upper_block_size,
        Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, Eigen::Dynamic>& H_p,
        Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, 1>& b_p,
        const typename DerivedMatrix::Scalar damping_factor,
        const typename DerivedMatrix::Scalar symmetry_violation_threshold, const bool jacobi_scaling);

/**
 * @brief Compute the square root \f$ S_p \f$, \f$ r_p \f$ of the Schur complement of the linear least squares problem
 * \f$ \min_{\delta x} \frac{1}{2} \| J \delta x + r \|^2 \f$, with \f$ J = [J_m \; J_k] \f$ and
 * \f$ \delta x = [\delta x_m^T \; \delta x_k^T]^T \f$, by QR decompositions which never form \f$ J^T J \f$, such that:
 * \f[
 *      \min_{\delta x_m} \frac{1}{2} \| J_m \delta x_m + J_k \delta x_k + r \|^2
 *          = \frac{1}{2} \| S_p \delta x_k + r_p \|^2 + c
 * \f]
 * for a constant \f$ c \f$. For the normal equations \f$ H \delta x = b \f$ with \f$ H = J^T J \f$ and
 * \f$ b = -J^T r \f$, this gives the `schur_complement` \f$ H_p = S_p^T S_p \f$ and \f$ b_p = -S_p^T r_p \f$, but with
 * an accuracy governed by the condition number \f$ \kappa(J) \f$ rather than \f$ \kappa(H) = \kappa(J)^2 \f$. The
 * round-off of forming \f$ H \f$ is of order \f$ \epsilon \| H \| \f$, which can exceed the information of weakly
 * constrained (e.g. gauge) directions, whereas here the error of their square root is of order \f$ \epsilon \| J \|
 * \f$.
 *
 * The column-pivoted QR decomposition \f$ J_m P_m = Q_m R_{mm} \f$ (\f$ P_m \f$ is the permutation matrix, \f$ Q_m \f$
 * is orthogonal, \f$ R_{mm} \f$ is upper triangular) of numerical rank \f$ p \f$ rotates the system:
 * \f[
 *      Q_m^T \begin{bmatrix} J_m P_m & J_k & r \end{bmatrix} =
 *      \begin{bmatrix} R_{mm} & \tilde{J}_{mk} & \tilde{r}_m \\ 0 & \tilde{J}_k & \tilde{r} \end{bmatrix}
 * \f]
 * whose first \f$ p \f$ rows are satisfied by \f$ \delta x_m \f$ for any \f$ \delta x_k \f$, leaving the marginal
 * \f$ \frac{1}{2} \| \tilde{J}_k \delta x_k + \tilde{r} \|^2 \f$ over the remaining rows. With the tolerance
 * \f$ \tau \f$ of `numerical_rank_tolerance` for the dimensions of \f$ J \f$, pivots at or below
 * \f$ \tau \max_i |R_{mm}(i, i)| \f$ are directions of \f$ \delta x_m \f$ unconstrained by \f$ J_m \f$, whose rows
 * stay in the marginal, so a rank-deficient \f$ J_m \f$ needs no damping. The column-pivoted QR decomposition
 * \f$ \tilde{J}_k P_p = Q_p R_p \f$ then compresses the marginal to at most the columns of \f$ J_k \f$ as rows:
 * \f[
 *      S_p = R_p P_p^T, \quad \begin{bmatrix} r_p \\ r_c \end{bmatrix} = Q_p^T \tilde{r}, \quad
 *      c = \frac{1}{2} \| r_c \|^2
 * \f]
 * Column pivoting bounds each row by its pivot, \f$ |R_p(i, j)| \leq |R_p(i, i)| \f$, so rows with a numerically zero
 * pivot, at or below \f$ \tau \max_i |R_p(i, i)| \f$ with the same tolerance, can be removed (e.g. with
 * `remove_rows_by_threshold`), discarding at most that magnitude. Both decompositions are no larger than \f$ J \f$, so
 * the tolerance for the dimensions of \f$ J \f$ bounds their round-off. The pivots behave like singular values, i.e.
 * the square roots of the eigenvalues of \f$ H_p \f$.
 *
 * The computation is in the scalar type of \f$ J \f$, so that e.g. `J.cast<long double>()` gives a reference of
 * higher precision, whose tolerance \f$ \tau \f$ then scales with the machine epsilon of that type.
 *
 * @tparam DerivedMatrix
 * @tparam DerivedVector
 * @param J Jacobian \f$ [J_m \; J_k] \f$
 * @param r residual vector
 * @param left_block_size columns of \f$ J_m \f$, or equivalently the column index of \f$ J_k \f$
 * @param S_p returned \f$ S_p \f$, with \f$ \min(\text{rows}(J) - p, \text{cols}(J_k)) \f$ rows and columns ordered as
 * \f$ J_k \f$
 * @param r_p returned \f$ r_p \f$, aligned with the rows of \f$ S_p \f$
 * @param pivots_mm returned absolute pivots \f$ |R_{mm}(i, i)| \f$ of \f$ J_m \f$, in non-increasing order
 * @param pivots_p returned absolute pivots \f$ |R_p(i, i)| \f$, in non-increasing order and aligned with the rows of
 * \f$ S_p \f$
 * @param rank_threshold_factor factor of the tolerance \f$ \tau \f$ of `numerical_rank_tolerance`, at or below which
 * relative to the maximum pivot a pivot of \f$ J_m \f$ is numerically zero. The default of 1 is the conventional
 * numerical rank tolerance, while a larger factor adds a margin, e.g. when the entries of \f$ J \f$ are inexact. A much
 * smaller factor risks treating a round-off pivot as genuine, which then projects a direction out of the marginal on
 * the basis of round-off. A much larger factor instead treats small but genuine pivots as zero, keeping rows in the
 * marginal that \f$ \delta x_m \f$ would satisfy, so that \f$ S_p^T S_p \f$ exceeds the exact Schur complement.
 * @return Eigen::Index numerical rank \f$ p \f$ of \f$ J_m \f$
 */
template<typename DerivedMatrix, typename DerivedVector>
    requires(std::is_same_v<typename DerivedMatrix::Scalar, typename DerivedVector::Scalar>)
Eigen::Index square_root_schur_complement(const Eigen::MatrixBase<DerivedMatrix>& J,
        const Eigen::MatrixBase<DerivedVector>& r, const int left_block_size,
        Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, Eigen::Dynamic>& S_p,
        Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, 1>& r_p,
        Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, 1>& pivots_mm,
        Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, 1>& pivots_p,
        const typename DerivedMatrix::Scalar rank_threshold_factor = static_cast<typename DerivedMatrix::Scalar>(1));

/**
 * @brief Compute the 3D skew-symmetric matrix \f$ [\mathbf{x}]_\times \f$ which when multiplied with a vector, is
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
 * Note also that \f$ \mathbf{x} \times \mathbf{y} = - \mathbf{y} \times \mathbf{x} = - [\mathbf{y}]_\times \mathbf{x} =
 * [\mathbf{y}]_\times^T \mathbf{x} = [\mathbf{x}]_\times \mathbf{y} \f$.
 *
 * @tparam Derived a 3x1 vector
 * @param v
 * @return Eigen::Matrix<typename Derived::Scalar, 3, 3>
 */
template<typename Derived>
    requires(Derived::RowsAtCompileTime == 3 && Derived::ColsAtCompileTime == 1)
constexpr Eigen::Matrix<typename Derived::Scalar, 3, 3> skew_symmetric_cross(const Eigen::MatrixBase<Derived>& v);

}

#include "mathbox/impl/matrix_operations.hpp"

#endif
