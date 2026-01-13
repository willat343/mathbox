#ifndef MATHBOX_IMPL_MATRIX_PROPERTIES_HPP
#define MATHBOX_IMPL_MATRIX_PROPERTIES_HPP

#include "mathbox/matrix_properties.hpp"

namespace math {

template<typename Derived>
inline bool has_positive_diagonals(const Eigen::MatrixBase<Derived>& m) {
    return m.diagonal().minCoeff() > 0.0;
}

template<typename Derived>
bool is_positive_definite(const typename Eigen::MatrixBase<Derived>& m) {
    if (!is_symmetric(m)) {
        return false;
    }
    typename Eigen::LLT<Derived> llt(m);  // compute the Cholesky decomposition
    return llt.info() == Eigen::Success;
}

template<typename Derived>
bool is_positive_semidefinite(const typename Eigen::MatrixBase<Derived>& m) {
    if (!is_symmetric(m)) {
        return false;
    }
    typename Eigen::SelfAdjointEigenSolver<Derived> solver(m);  // compute the Cholesky decomposition
    if (solver.info() != Eigen::Success) {
        return false;
    }
    return (solver.eigenvalues().array() >= static_cast<typename Derived::Scalar>(0)).all();
}

template<typename Derived>
inline bool is_skew_symmetric(const Eigen::DenseBase<Derived>& m) {
    return m.isApprox(-m.transpose());
}

template<typename Derived>
bool is_symmetric(const Eigen::DenseBase<Derived>& m, const typename Derived::Scalar precision) {
    for (Eigen::Index r = 0; r < m.rows() - 1; ++r) {
        for (Eigen::Index c = r + 1; c < m.cols(); ++c) {
            if (std::abs(m(r, c) - m(c, r)) > precision) {
                return false;
            }
        }
    }
    return true;
}

template<typename Derived>
bool is_upper_triangular(const Eigen::DenseBase<Derived>& m, const typename Derived::Scalar precision) {
    for (Eigen::Index r = 0; r < m.rows() - 1; ++r) {
        for (Eigen::Index c = 0; c < r; ++c) {
            if (std::abs(m(r, c)) > precision) {
                return false;
            }
        }
    }
    return true;
}

template<typename Derived>
int num_zero_columns(const Eigen::MatrixBase<Derived>& m, const typename Derived::Scalar precision) {
    std::size_t num_zero_columns_{0};
    for (int i = 0; i < m.cols(); ++i) {
        num_zero_columns_ += (m.col(i).isZero(precision) ? 1 : 0);
    }
    return num_zero_columns_;
}

template<typename Derived>
int num_zero_rows(const Eigen::MatrixBase<Derived>& m, const typename Derived::Scalar precision) {
    std::size_t num_zero_rows_{0};
    for (int i = 0; i < m.rows(); ++i) {
        num_zero_rows_ += (m.row(i).isZero(precision) ? 1 : 0);
    }
    return num_zero_rows_;
}

}

#endif
