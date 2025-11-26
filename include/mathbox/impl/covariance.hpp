#ifndef MATHBOX_IMPL_COVARIANCE_HPP
#define MATHBOX_IMPL_COVARIANCE_HPP

#include <cppbox/exceptions.hpp>

#include "mathbox/covariance.hpp"
#include "mathbox/matrix_properties.hpp"

namespace math {

template<typename Scalar_, int Size_>
template<typename Derived>
PositiveSemiDefiniteMatrix<Scalar_, Size_>::PositiveSemiDefiniteMatrix(const Eigen::MatrixBase<Derived>& matrix,
        const bool skip_checks)
    requires(SizeAtCompileTime == 1 || Derived::ColsAtCompileTime != 1)
    : Base(matrix) {
    throw_if(!skip_checks && !math::is_positive_semidefinite(matrix.eval()),
            "Cannnot initialise PositiveSemiDefiniteMatrix with matrix that isn't positive semi-definite.");
}

template<typename Scalar_, int Size_>
template<typename Derived>
PositiveSemiDefiniteMatrix<Scalar_, Size_>::PositiveSemiDefiniteMatrix(const Eigen::MatrixBase<Derived>& diagonal,
        const bool skip_checks)
    requires(SizeAtCompileTime != 1 && Derived::ColsAtCompileTime == 1)
    : Base(Matrix(diagonal.asDiagonal())) {
    throw_if(!skip_checks && diagonal.minCoeff() < static_cast<Scalar>(0),
            "Cannot initialise PositiveSemiDefiniteMatrix with any negative diagonal elements.");
}

template<typename Scalar_, int Size_>
PositiveSemiDefiniteMatrix<Scalar_, Size_>::PositiveSemiDefiniteMatrix(const Scalar diagonal_element, const int size,
        const bool skip_checks)
    : Base(Matrix(Vector::Constant(size, diagonal_element).asDiagonal())) {
    if (!skip_checks) {
        throw_if(diagonal_element < static_cast<Scalar>(0),
                "Cannot initialise PositiveSemiDefiniteMatrix with negative single diagonal element.");
        if constexpr (SizeAtCompileTime == Eigen::Dynamic) {
            throw_if(size == Eigen::Dynamic,
                    "Cannot initialise dynamic-size PositiveSemiDefiniteMatrix with single diagonal element without "
                    "explicitly setting size (size was Eigen::Dynamic).");
        } else {
            throw_if(size != SizeAtCompileTime,
                    "Cannot initialise fixed-size PositiveSemiDefiniteMatrix with single diagonal "
                    "element and a size that does not match the fixed size.");
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

template<typename Derived>
Eigen::Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, Derived::RowsAtCompileTime> covariance_from_sigmas(
        const Eigen::MatrixBase<Derived>& sigmas) {
    if constexpr (Derived::ColsAtCompileTime == Eigen::Dynamic) {
        throw_if(sigmas.cols() != 1, "Sigmas was not a column vector.");
    } else {
        static_assert(Derived::ColsAtCompileTime == 1, "Sigmas was not a column vector.");
    }
    return covariance_from_sigmas<Derived::RowsAtCompileTime, typename Derived::Scalar>(sigmas);
}

template<int Rows, typename Scalar>
inline Eigen::Matrix<Scalar, Rows, Rows> covariance_from_sigmas(
        const Eigen::Ref<const Eigen::Matrix<Scalar, Rows, 1>>& sigmas) {
    return covariance_from_variances<Rows, Scalar>(sigmas.cwiseProduct(sigmas));
}

template<int Rows, typename Scalar>
inline Eigen::Matrix<Scalar, Rows, Rows> covariance_from_sigma(const Scalar sigma) {
    static_assert(Rows > 0, "number of compile-time rows must be > 0 when constructing from a single sigma value");
    return Eigen::Matrix<Scalar, Rows, 1>::Constant(covariance_from_sigma(sigma)).asDiagonal();
}

template<typename Scalar>
Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> covariance_from_sigma(const Scalar sigma, const int rows) {
    throw_if(rows < 1, "covariance_from_sigma expected rows to be a positive integer");
    return Eigen::Matrix<Scalar, Eigen::Dynamic, 1>::Constant(rows, covariance_from_sigma(sigma)).asDiagonal();
}

template<typename Scalar>
inline Scalar covariance_from_sigma(const Scalar sigma) {
    return sigma * sigma;
}

template<typename Derived>
Eigen::Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, Derived::RowsAtCompileTime>
covariance_from_variances(const Eigen::MatrixBase<Derived>& variances) {
    if constexpr (Derived::ColsAtCompileTime == Eigen::Dynamic) {
        throw_if(variances.cols() != 1, "Variances was not a column vector.");
    } else {
        static_assert(Derived::ColsAtCompileTime == 1, "Variances was not a column vector.");
    }
    return covariance_from_variances<Derived::RowsAtCompileTime, typename Derived::Scalar>(variances);
}

template<int Rows, typename Scalar>
inline Eigen::Matrix<Scalar, Rows, Rows> covariance_from_variances(
        const Eigen::Ref<const Eigen::Matrix<Scalar, Rows, 1>>& variances) {
    return variances.asDiagonal();
}

template<int Rows, typename Scalar>
inline Eigen::Matrix<Scalar, Rows, Rows> covariance_from_variance(const Scalar variance) {
    static_assert(Rows > 0, "number of compile-time rows must be > 0 when constructing from a single variance value");
    return Eigen::Matrix<Scalar, Rows, 1>::Constant(covariance_from_variance(variance)).asDiagonal();
}

template<typename Scalar>
inline Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> covariance_from_variance(const Scalar variance,
        const int rows) {
    throw_if(rows < 1, "covariance_from_variance expected rows to be a positive integer");
    return Eigen::Matrix<Scalar, Eigen::Dynamic, 1>::Constant(rows, covariance_from_variance(variance)).asDiagonal();
}

template<typename Scalar>
inline Scalar covariance_from_variance(const Scalar variance) {
    return variance;
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

#if !MATHBOX_HEADER_ONLY
namespace math {

extern template class PositiveSemiDefiniteMatrix<double, Eigen::Dynamic>;
extern template class PositiveSemiDefiniteMatrix<double, 1>;
extern template class PositiveSemiDefiniteMatrix<double, 2>;
extern template class PositiveSemiDefiniteMatrix<double, 3>;
extern template class PositiveSemiDefiniteMatrix<double, 4>;
extern template class PositiveSemiDefiniteMatrix<double, 5>;
extern template class PositiveSemiDefiniteMatrix<double, 6>;
extern template class PositiveSemiDefiniteMatrix<double, 7>;
extern template class PositiveSemiDefiniteMatrix<double, 8>;
extern template class PositiveSemiDefiniteMatrix<double, 9>;

extern template class CovarianceDensity<double, Eigen::Dynamic>;
extern template class CovarianceDensity<double, 1>;
extern template class CovarianceDensity<double, 2>;
extern template class CovarianceDensity<double, 3>;
extern template class CovarianceDensity<double, 4>;
extern template class CovarianceDensity<double, 5>;
extern template class CovarianceDensity<double, 6>;
extern template class CovarianceDensity<double, 7>;
extern template class CovarianceDensity<double, 8>;
extern template class CovarianceDensity<double, 9>;

}
#endif

#endif
