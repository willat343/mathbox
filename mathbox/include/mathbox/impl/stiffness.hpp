#ifndef MATHBOX_IMPL_STIFFNESS_HPP
#define MATHBOX_IMPL_STIFFNESS_HPP

#include <cmath>

#include "mathbox/decompose.hpp"
#include "mathbox/stiffness.hpp"

namespace math {

template<typename Derived>
inline Derived stiffness_from_covariance(const Eigen::MatrixBase<Derived>& covariance) {
    if constexpr (Derived::RowsAtCompileTime == Eigen::Dynamic || Derived::ColsAtCompileTime == Eigen::Dynamic) {
        if (covariance.rows() != covariance.cols()) {
            throw std::runtime_error("Covariance is not square.");
        }
    } else {
        static_assert(Derived::RowsAtCompileTime == Derived::ColsAtCompileTime, "Covariance is not square.");
    }
    return UTU(covariance.inverse().eval(), math::LLTDecompositionMethod::CHOLESKY);
}

template<int Rows, typename Scalar>
Eigen::Matrix<Scalar, Rows, Rows> stiffness_from_covariance(
        const Eigen::Ref<const Eigen::Matrix<Scalar, Rows, Rows>>& covariance) {
    return stiffness_from_covariance<Eigen::Matrix<Scalar, Rows, Rows>>(covariance.eval());
}

template<typename Derived>
inline Derived stiffness_from_information(const Eigen::MatrixBase<Derived>& information) {
    if constexpr (Derived::RowsAtCompileTime == Eigen::Dynamic || Derived::ColsAtCompileTime == Eigen::Dynamic) {
        if (information.rows() != information.cols()) {
            throw std::runtime_error("Information is not square.");
        }
    } else {
        static_assert(Derived::RowsAtCompileTime == Derived::ColsAtCompileTime, "Information is not square.");
    }
    return UTU(information, math::LLTDecompositionMethod::CHOLESKY);
}

template<int Rows, typename Scalar>
Eigen::Matrix<Scalar, Rows, Rows> stiffness_from_information(
        const Eigen::Ref<const Eigen::Matrix<Scalar, Rows, Rows>>& information) {
    return stiffness_from_information<Eigen::Matrix<Scalar, Rows, Rows>>(information.eval());
}

template<typename Derived>
Eigen::Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, Derived::RowsAtCompileTime> stiffness_from_sigmas(
        const Eigen::MatrixBase<Derived>& sigmas) {
    if constexpr (Derived::ColsAtCompileTime == Eigen::Dynamic) {
        if (sigmas.cols() != 1) {
            throw std::runtime_error("Sigmas was not a column vector.");
        }
    } else {
        static_assert(Derived::ColsAtCompileTime == 1, "Sigmas was not a column vector.");
    }
    return stiffness_from_sigmas<Derived::RowsAtCompileTime, typename Derived::Scalar>(sigmas);
}

template<int Rows, typename Scalar>
inline Eigen::Matrix<Scalar, Rows, Rows> stiffness_from_sigmas(
        const Eigen::Ref<const Eigen::Matrix<Scalar, Rows, 1>>& sigmas) {
    return sigmas.cwiseInverse().asDiagonal();
}

template<int Rows, typename Scalar>
inline Eigen::Matrix<Scalar, Rows, Rows> stiffness_from_sigma(const Scalar sigma) {
    static_assert(Rows > 0, "number of compile-time rows must be > 0 when constructing from a single sigma value");
    return Eigen::Matrix<Scalar, Rows, 1>::Constant(stiffness_from_sigma(sigma)).asDiagonal();
}

template<typename Scalar>
Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> stiffness_from_sigma(const Scalar sigma, const int rows) {
    return Eigen::Matrix<Scalar, Eigen::Dynamic, 1>::Constant(rows, 1, stiffness_from_sigma(sigma)).asDiagonal();
}

template<typename Scalar>
inline Scalar stiffness_from_sigma(const Scalar sigma) {
    return static_cast<Scalar>(1) / sigma;
}

template<typename Derived>
Eigen::Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, Derived::RowsAtCompileTime>
stiffness_from_variances(const Eigen::MatrixBase<Derived>& variances) {
    if constexpr (Derived::ColsAtCompileTime == Eigen::Dynamic) {
        if (variances.cols() != 1) {
            throw std::runtime_error("Variances was not a column vector.");
        }
    } else {
        static_assert(Derived::ColsAtCompileTime == 1, "Variances was not a column vector.");
    }
    return stiffness_from_variances<Derived::RowsAtCompileTime, typename Derived::Scalar>(variances);
}

template<int Rows, typename Scalar>
inline Eigen::Matrix<Scalar, Rows, Rows> stiffness_from_variances(
        const Eigen::Ref<const Eigen::Matrix<Scalar, Rows, 1>>& variances) {
    return stiffness_from_sigmas<Rows, Scalar>(variances.cwiseSqrt());
}

template<int Rows, typename Scalar>
inline Eigen::Matrix<Scalar, Rows, Rows> stiffness_from_variance(const Scalar variance) {
    static_assert(Rows > 0, "number of compile-time rows must be > 0 when constructing from a single variance value");
    return Eigen::Matrix<Scalar, Rows, 1>::Constant(stiffness_from_variance(variance)).asDiagonal();
}

template<typename Scalar>
inline Scalar stiffness_from_variance(const Scalar variance) {
    return static_cast<Scalar>(1) / std::sqrt(variance);
}

}

#endif
