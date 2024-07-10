#ifndef MATHBOX_IMPL_MATRIX_PROPERTIES_HPP
#define MATHBOX_IMPL_MATRIX_PROPERTIES_HPP

#include "mathbox/matrix_properties.hpp"

namespace math {

template<typename Derived>
bool has_positive_diagonals(const Eigen::MatrixBase<Derived>& m) {
    return m.diagonal().minCoeff() > 0.0;
}

template<typename Derived>
bool is_positive_definite(const typename Eigen::EigenBase<Derived>& m) {
    if (!is_symmetric(m)) {
        return false;
    }
    typename Eigen::LLT<Derived> llt(m);  // compute the Cholesky decomposition
    return llt.info() == Eigen::Success;
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

}

#endif
