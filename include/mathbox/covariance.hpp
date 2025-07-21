#ifndef MATHBOX_COVARIANCE_HPP
#define MATHBOX_COVARIANCE_HPP

#include <Eigen/Core>
#include <type_traits>

namespace math {

/**
 * @brief PositiveSemiDefiniteMatrix class
 *
 * @tparam Scalar_
 * @tparam Size_ positive integer or `Eigen::Dynamic`
 */
template<typename Scalar_, int Size_>
class PositiveSemiDefiniteMatrix : public Eigen::Matrix<Scalar_, Size_, Size_> {
public:
    using Scalar = Scalar_;
    static constexpr int SizeAtCompileTime = Size_;
    using Base = Eigen::Matrix<Scalar_, Size_, Size_>;
    using Matrix = Base;
    using Vector = Eigen::Matrix<Scalar_, Size_, 1>;

    /**
     * @brief Construct a Positive Semi-Definite Matrix from a matrix.
     *
     * @tparam Derived
     * @param matrix
     * @param skip_checks
     *
     * @tparam Derived
     */
    template<typename Derived>
    explicit PositiveSemiDefiniteMatrix(const Eigen::MatrixBase<Derived>& matrix, const bool skip_checks = false)
        requires(SizeAtCompileTime == 1 || Derived::ColsAtCompileTime != 1);

    /**
     * @brief Construct a Positive Semi-Definite Matrix from a diagonal (column) vector.
     *
     * All off-diagonal elements are initialised to zero.
     *
     * @tparam Derived
     * @param diagonal diagonal elements (e.g. variances for a Covariance)
     * @param skip_checks
     */
    template<typename Derived>
    explicit PositiveSemiDefiniteMatrix(const Eigen::MatrixBase<Derived>& diagonal, const bool skip_checks = false)
        requires(SizeAtCompileTime != 1 && Derived::ColsAtCompileTime == 1);

    /**
     * @brief Construct a dynamic-size (or fixed-size) Positive Semi-Definite Matrix object with all elements along the
     * diagonal equal to `diagonal_element`.
     *
     * All off-diagonal elements are initialised to zero.
     *
     * For fixed-size matrices, size is checked against `SizeAtCompileTime`.
     *
     * @param diagonal_element diagonal element (e.g. variance for a Covariance)
     * @param size size of matrix
     * @param skip_checks
     */
    explicit PositiveSemiDefiniteMatrix(const Scalar diagonal_element, const int size, const bool skip_checks = false);

    /**
     * @brief Construct a fixed-size Positive Semi-Definite Matrix with all elements along the diagonal equal to
     * `diagonal_element`.
     *
     * All off-diagonal elements are initialised to zero.
     *
     * @param diagonal_element diagonal element (e.g. variance for a Covariance)
     * @param skip_checks
     */
    explicit PositiveSemiDefiniteMatrix(const Scalar diagonal_element, const bool skip_checks = false);

    /**
     * @brief Assign from a positive semi-definite matrix or a diagonal (column) vector.
     *
     * @tparam Derived
     */
    template<typename Derived>
    PositiveSemiDefiniteMatrix& operator=(const Eigen::MatrixBase<Derived>& rhs);
};

template<typename Scalar, int Size>
using Covariance = PositiveSemiDefiniteMatrix<Scalar, Size>;
using CovarianceXd = Covariance<double, Eigen::Dynamic>;
using Covariance1d = Covariance<double, 1>;
using Covariance2d = Covariance<double, 2>;
using Covariance3d = Covariance<double, 3>;
using Covariance4d = Covariance<double, 4>;
using Covariance5d = Covariance<double, 5>;
using Covariance6d = Covariance<double, 6>;
using Covariance7d = Covariance<double, 7>;
using Covariance8d = Covariance<double, 8>;
using Covariance9d = Covariance<double, 9>;
template<typename Scalar>
using CovarianceX = Covariance<Scalar, Eigen::Dynamic>;
template<typename Scalar>
using Covariance1 = Covariance<Scalar, 1>;
template<typename Scalar>
using Covariance2 = Covariance<Scalar, 2>;
template<typename Scalar>
using Covariance3 = Covariance<Scalar, 3>;
template<typename Scalar>
using Covariance4 = Covariance<Scalar, 4>;
template<typename Scalar>
using Covariance5 = Covariance<Scalar, 5>;
template<typename Scalar>
using Covariance6 = Covariance<Scalar, 6>;
template<typename Scalar>
using Covariance7 = Covariance<Scalar, 7>;
template<typename Scalar>
using Covariance8 = Covariance<Scalar, 8>;
template<typename Scalar>
using Covariance9 = Covariance<Scalar, 9>;
template<int Size>
using Covarianced = Covariance<double, Size>;

/**
 * @brief Covariance density class, which can be considered as a covariance per unit (e.g. temporal, spatial, spectral,
 * or other).
 *
 * This class inherits all the runtime requirements of a Positive Semi-Definite Matrix, but additionally provides
 * functionally to compute covariance functions.
 *
 * @tparam Scalar_
 * @tparam Size_ positive integer or `Eigen::Dynamic`
 */
template<typename Scalar_, int Size_>
class CovarianceDensity : public PositiveSemiDefiniteMatrix<Scalar_, Size_> {
public:
    using Scalar = Scalar_;
    static constexpr int SizeAtCompileTime = Size_;
    using Base = Covariance<Scalar_, Size_>;
    using Matrix = Base;
    using Vector = Eigen::Matrix<Scalar_, Size_, 1>;

    template<typename Derived>
    explicit CovarianceDensity(const Eigen::MatrixBase<Derived>& covariance_density_or_variances_density,
            const bool skip_checks = false);

    explicit CovarianceDensity(const Scalar variance_density, const int size, const bool skip_checks = false);

    explicit CovarianceDensity(const Scalar variance_density, const bool skip_checks = false);

    /**
     * @brief Compute a covariance from the covariance density by multiplying by `span`, which can be considered a
     * duration in a temporal context, or length in a spatial context.
     *
     * Note that checks are skipped in the creation of covariances, but if this class was constructed with `skip_checks`
     * disabled then the produced covariances are guaranteed to be valid.
     *
     * @param span
     * @return Covariance<Scalar, SizeAtCompileTime>
     */
    Covariance<Scalar, SizeAtCompileTime> covariance(const Scalar span) const;

    /**
     * @brief Assign from a covariance matrix or variances vector.
     *
     * @tparam Derived
     */
    template<typename Derived>
    CovarianceDensity& operator=(const Eigen::MatrixBase<Derived>& rhs);
};

using CovarianceDensityXd = CovarianceDensity<double, Eigen::Dynamic>;
using CovarianceDensity1d = CovarianceDensity<double, 1>;
using CovarianceDensity2d = CovarianceDensity<double, 2>;
using CovarianceDensity3d = CovarianceDensity<double, 3>;
using CovarianceDensity4d = CovarianceDensity<double, 4>;
using CovarianceDensity5d = CovarianceDensity<double, 5>;
using CovarianceDensity6d = CovarianceDensity<double, 6>;
using CovarianceDensity7d = CovarianceDensity<double, 7>;
using CovarianceDensity8d = CovarianceDensity<double, 8>;
using CovarianceDensity9d = CovarianceDensity<double, 9>;
template<typename Scalar>
using CovarianceDensityX = CovarianceDensity<Scalar, Eigen::Dynamic>;
template<typename Scalar>
using CovarianceDensity1 = CovarianceDensity<Scalar, 1>;
template<typename Scalar>
using CovarianceDensity2 = CovarianceDensity<Scalar, 2>;
template<typename Scalar>
using CovarianceDensity3 = CovarianceDensity<Scalar, 3>;
template<typename Scalar>
using CovarianceDensity4 = CovarianceDensity<Scalar, 4>;
template<typename Scalar>
using CovarianceDensity5 = CovarianceDensity<Scalar, 5>;
template<typename Scalar>
using CovarianceDensity6 = CovarianceDensity<Scalar, 6>;
template<typename Scalar>
using CovarianceDensity7 = CovarianceDensity<Scalar, 7>;
template<typename Scalar>
using CovarianceDensity8 = CovarianceDensity<Scalar, 8>;
template<typename Scalar>
using CovarianceDensity9 = CovarianceDensity<Scalar, 9>;
template<int Size>
using CovarianceDensityd = CovarianceDensity<double, Size>;

}

#include "mathbox/impl/covariance.hpp"

#endif
