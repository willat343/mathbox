#ifndef MATHBOX_IMPL_DECOMPOSE_HPP
#define MATHBOX_IMPL_DECOMPOSE_HPP

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>

#include "mathbox/decompose.hpp"

namespace math {

inline void check_llt_computation_info(const Eigen::ComputationInfo info) {
    if (info != Eigen::ComputationInfo::Success) {
        if (info == Eigen::ComputationInfo::NumericalIssue) {
            throw std::runtime_error("NumericalIssue encountered while performing LL^* decomposition of matrix.");
        } else if (info == Eigen::ComputationInfo::NoConvergence) {
            throw std::runtime_error("NoConvergence encountered while performing LL^* decomposition of matrix.");
        } else if (info == Eigen::InvalidInput) {
            throw std::runtime_error("InvalidInput encountered while performing LL^* decomposition of matrix.");
        } else {
            throw std::runtime_error("Unknown error encountered while performing LL^* decomposition of matrix.");
        }
    }
}

template<typename Derived>
Derived LLT(const Eigen::MatrixBase<Derived>& matrix, const LLTDecompositionMethod method) {
    Derived L;
    if (method == LLTDecompositionMethod::CHOLESKY) {
        // Decomposition is LL^T
        Eigen::LLT<Derived> decomposer(matrix);
        check_llt_computation_info(decomposer.info());
        L = decomposer.matrixL();
    } else if (method == LLTDecompositionMethod::EIGEN) {
        // Decomposition is QDQ^T
        Eigen::SelfAdjointEigenSolver<Derived> decomposer(matrix);
        check_llt_computation_info(decomposer.info());
        L = decomposer.eigenvectors() * decomposer.eigenvalues().cwiseSqrt().asDiagonal();
    } else if (method == LLTDecompositionMethod::ROBUST_CHOLESKY) {
        // Decomposition is P^TLDL^TP
        Eigen::LDLT<Derived> decomposer(matrix);
        check_llt_computation_info(decomposer.info());
        L = decomposer.vectorD().cwiseSqrt().asDiagonal();
        L = decomposer.matrixL() * L;
        L = decomposer.transpositionsP().transpose() * L;
    } else {
        throw std::runtime_error("DecompositionMethod " + std::to_string(method) + " not recognised.");
    }
    return L;
}

template<typename Derived>
inline Derived UTU(const Eigen::MatrixBase<Derived>& matrix, const LLTDecompositionMethod method) {
    return LLT(matrix, method).transpose();
}

}

#endif
