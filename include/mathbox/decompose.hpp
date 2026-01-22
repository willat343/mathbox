#ifndef MATHBOX_DECOMPOSE_HPP
#define MATHBOX_DECOMPOSE_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>

namespace math {

/**
 * @brief LLT or UTU decomposition method for a matrix.
 *
 */
enum LLTDecompositionMethod {
    CHOLESKY,       /**< Cholesky decomposition. Usually faster than Eigen decomposition, but valid only for positive
                        definite matrices, not positive semi-definite matrices. Returns the lower triangular matrix. */
    EIGEN,          /**< Eigen decomposition. Only valid for unitary (orthogonal), Hermitian (symmetric) or
                        skew-Hermitian (skew-symmetric) matrices (e.g. covariance/information matrices). Note that
                        neither the element nor their ordering is in general the same as any other method. Note also
                        that eigenvectors and eigenvalues are obtained ordered, which determines the ordering of L. */
    ROBUST_CHOLESKY /**< Robust Cholesky decomposition. Valid for positive semi-definite and negative semi-definite
                        matrices. Note that neither the element nor their ordering is in general the same as any other
                        method. */
};

/**
 * @brief List of all decomposition methods.
 *
 */
inline const std::vector<LLTDecompositionMethod> llt_decomposition_methods = {{CHOLESKY, EIGEN, ROBUST_CHOLESKY}};

/**
 * @brief Decompose a hermitian, matrix (of valid form, see `LLTDecompositionMethod` documentation) into form
 * \f$\mathbf{L}\mathbf{L}^*\f$ (\f$\mathbf{L}\mathbf{L}^T\f$ for real matrices), and return \f$\mathbf{L}\f$.
 *
 * Note that different LLT decomposition methods will yield different \f$L\f$ matrices.
 *
 * @tparam Derived
 * @param matrix matrix to decompose
 * @param method the LLT decomposition method
 * @return Derived \f$L\f$ matrix
 */
template<typename Derived>
Derived LLT(const Eigen::MatrixBase<Derived>& matrix, const LLTDecompositionMethod method = ROBUST_CHOLESKY);

/**
 * @brief Decompose a hermitian, matrix (of valid form, see `math::LLTDecompositionMethod`) into form
 * \f$\mathbf{U}^*\mathbf{U}\f$ (\f$\mathbf{U}^T\mathbf{U}\f$ for real matrices), and return \f$\mathbf{U}\f$.
 *
 * Note that different LLT decomposition methods will yield different \f$L\f$ matrices.
 *
 * @tparam Derived
 * @param matrix matrix to decompose
 * @param method the LLT decomposition method
 * @return Derived \f$U\f$ matrix
 */
template<typename Derived>
Derived UTU(const Eigen::MatrixBase<Derived>& matrix, const LLTDecompositionMethod method = ROBUST_CHOLESKY);

}

#include "mathbox/impl/decompose.hpp"

#endif
