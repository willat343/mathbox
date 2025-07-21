#ifndef MATHBOX_IMPL_COVARIANCE_HPP
#define MATHBOX_IMPL_COVARIANCE_HPP

#include "mathbox/covariance.hpp"
#include "mathbox/matrix_properties.hpp"

namespace math {

template<typename Scalar_, int Size_>
template<typename Derived>
PositiveSemiDefiniteMatrix<Scalar_, Size_>::PositiveSemiDefiniteMatrix(const Eigen::MatrixBase<Derived>& matrix,
        const bool skip_checks)
    requires(SizeAtCompileTime == 1 || Derived::ColsAtCompileTime != 1)
    : Base(matrix) {
    if (!skip_checks && !math::is_positive_semidefinite(matrix.eval())) {
        throw std::runtime_error(
                "Cannnot initialise PositiveSemiDefiniteMatrix with matrix that isn't positive semi-definite.");
    }
}

template<typename Scalar_, int Size_>
template<typename Derived>
PositiveSemiDefiniteMatrix<Scalar_, Size_>::PositiveSemiDefiniteMatrix(const Eigen::MatrixBase<Derived>& diagonal,
        const bool skip_checks)
    requires(SizeAtCompileTime != 1 && Derived::ColsAtCompileTime == 1)
    : Base(Matrix(diagonal.asDiagonal())) {
    if (!skip_checks && diagonal.minCoeff() < static_cast<Scalar>(0)) {
        throw std::runtime_error("Cannot initialise PositiveSemiDefiniteMatrix with any negative diagonal elements.");
    }
}

template<typename Scalar_, int Size_>
PositiveSemiDefiniteMatrix<Scalar_, Size_>::PositiveSemiDefiniteMatrix(const Scalar diagonal_element, const int size,
        const bool skip_checks)
    : Base(Matrix(Vector::Constant(size, diagonal_element).asDiagonal())) {
    if (!skip_checks) {
        if (diagonal_element < static_cast<Scalar>(0)) {
            throw std::runtime_error(
                    "Cannot initialise PositiveSemiDefiniteMatrix with negative single diagonal element.");
        }
        if constexpr (SizeAtCompileTime == Eigen::Dynamic) {
            if (size == Eigen::Dynamic) {
                throw std::runtime_error("Cannot initialise dynamic-size PositiveSemiDefiniteMatrix with single "
                                         "diagonal element without explicitly setting size (size was Eigen::Dynamic).");
            }
        } else {
            if (size != SizeAtCompileTime) {
                throw std::runtime_error("Cannot initialise fixed-size PositiveSemiDefiniteMatrix with single diagonal "
                                         "element and a size that does not match the fixed size.");
            }
        }
    }
}

template<typename Scalar_, int Size_>
PositiveSemiDefiniteMatrix<Scalar_, Size_>::PositiveSemiDefiniteMatrix(const Scalar diagonal_element,
        const bool skip_checks)
    : PositiveSemiDefiniteMatrix(diagonal_element, SizeAtCompileTime, skip_checks) {}

template<typename Scalar_, int Size_>
template<typename Derived>
PositiveSemiDefiniteMatrix<Scalar_, Size_>& PositiveSemiDefiniteMatrix<Scalar_, Size_>::operator=(
        const Eigen::MatrixBase<Derived>& rhs) {
    *this = PositiveSemiDefiniteMatrix{rhs, false};
    return *this;
}

template<typename Scalar_, int Size_>
template<typename Derived>
CovarianceDensity<Scalar_, Size_>::CovarianceDensity(
        const Eigen::MatrixBase<Derived>& covariance_density_or_variances_density, const bool skip_checks)
    : Base(covariance_density_or_variances_density, skip_checks) {}

template<typename Scalar_, int Size_>
CovarianceDensity<Scalar_, Size_>::CovarianceDensity(const Scalar variance, const int size, const bool skip_checks)
    : Base(variance, size, skip_checks) {}

template<typename Scalar_, int Size_>
CovarianceDensity<Scalar_, Size_>::CovarianceDensity(const Scalar variance, const bool skip_checks)
    : Base(variance, skip_checks) {}

template<typename Scalar_, int Size_>
inline auto CovarianceDensity<Scalar_, Size_>::covariance(const Scalar span) const
        -> Covariance<Scalar, SizeAtCompileTime> {
    return Covariance<Scalar, SizeAtCompileTime>(*this * span, true);
}

template<typename Scalar_, int Size_>
template<typename Derived>
CovarianceDensity<Scalar_, Size_>& CovarianceDensity<Scalar_, Size_>::operator=(const Eigen::MatrixBase<Derived>& rhs) {
    *this = CovarianceDensity{rhs, false};
    return *this;
}

}

#endif
