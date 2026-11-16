#ifndef MATHBOX_IMPL_DIFFERENTIATOR_HPP
#define MATHBOX_IMPL_DIFFERENTIATOR_HPP

#include <Eigen/QR>
#include <algorithm>
#include <cmath>
#include <cppbox/exceptions.hpp>
#include <cppbox/time.hpp>
#include <functional>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

#include "mathbox/combinatorics.hpp"
#include "mathbox/differentiator.hpp"

namespace math {

template<std::floating_point Scalar_>
inline unsigned int LocalPolynomialDifferentiator<Scalar_>::derivative_order() const {
    return derivative_order_;
}

template<std::floating_point Scalar_>
inline std::size_t LocalPolynomialDifferentiator<Scalar_>::evaluation_index() const {
    return evaluation_index_;
}

template<std::floating_point Scalar_>
inline unsigned int LocalPolynomialDifferentiator<Scalar_>::polynomial_degree() const {
    return polynomial_degree_;
}

template<std::floating_point Scalar_>
inline std::size_t LocalPolynomialDifferentiator<Scalar_>::window_size() const {
    return window_size_;
}

template<std::floating_point Scalar_>
LocalPolynomialDifferentiator<Scalar_>::LocalPolynomialDifferentiator(const std::size_t window_size_,
        const unsigned int polynomial_degree_, const unsigned int derivative_order_,
        const std::size_t evaluation_index_)
    : window_size_(window_size_),
      polynomial_degree_(polynomial_degree_),
      derivative_order_(derivative_order_),
      evaluation_index_(evaluation_index_) {
    // Check that the samples determine the polynomial (N > p), that its derivative is not zero by construction
    // (d <= p) and that the evaluated sample is in the window
    throw_if(window_size() <= polynomial_degree(), "Window size must be greater than the polynomial degree.",
            std::invalid_argument);
    throw_if(derivative_order() > polynomial_degree(), "Derivative order must not exceed the polynomial degree.",
            std::invalid_argument);
    throw_if(evaluation_index() >= window_size(), "Evaluation index must be within the window.", std::invalid_argument);
}

template<std::floating_point Scalar_>
inline std::size_t LocalPolynomialDifferentiator<Scalar_>::central_index(const std::size_t window_size) {
    throw_if(window_size % 2 == 0, "Window size must be odd for its centre to be a sample.", std::invalid_argument);
    return (window_size - 1) / 2;
}

template<std::floating_point Scalar_>
template<std::ranges::forward_range Range>
inline void LocalPolynomialDifferentiator<Scalar_>::check_range_size(const Range& range) const {
    throw_if(std::cmp_not_equal(std::ranges::distance(range), window_size()),
            "Range size (" + std::to_string(std::ranges::distance(range)) + ") must equal the window size (" +
                    std::to_string(window_size()) + ").",
            std::invalid_argument);
}

template<std::floating_point Scalar_>
auto LocalPolynomialDifferentiator<Scalar_>::compute_weights(
        const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& offsets) const -> Eigen::Matrix<Scalar, Eigen::Dynamic, 1> {
    // Check that the offsets are finite and strictly increasing, so distinct, where the window size, polynomial degree
    // and derivative order were checked on construction
    throw_if(!offsets.allFinite() ||
                     std::adjacent_find(offsets.begin(), offsets.end(), std::greater_equal<Scalar>()) != offsets.end(),
            "Sample points must be finite and strictly increasing.", std::invalid_argument);

    // Scale the offsets to xi_i = delta_i / h in [-1, 1] with h = max_i |delta_i|, so that their powers are in [-1, 1]
    // for any unit of x (the accuracy barely depends on it, since Householder QR is insensitive to the scale of the
    // columns), except for a single sample at the evaluation point, where h = 1 since any scale gives xi_0 = 0
    const Scalar max_absolute_offset = offsets.cwiseAbs().maxCoeff();
    const Scalar scale = max_absolute_offset > static_cast<Scalar>(0) ? max_absolute_offset : static_cast<Scalar>(1);
    // Vandermonde matrix of size N x (p + 1), A_ij = xi_i^j, built column by column, A_ij = A_i(j-1) xi_i
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> vandermonde(offsets.size(), polynomial_degree() + 1);
    vandermonde.col(0).setOnes();
    for (unsigned int j = 1; j <= polynomial_degree(); ++j) {
        vandermonde.col(j) = vandermonde.col(j - 1).cwiseProduct(offsets / scale);
    }

    // Fit the polynomial y(xi) = sum_j c_j xi^j to the samples y in the least squares sense with the column-pivoted QR
    // decomposition A P = Q R, c = P R_1^-1 Q_1^T y, where the leading block R_1 of R is invertible since the N > p
    // offsets are distinct
    const Eigen::ColPivHouseholderQR<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>> vandermonde_qr(vandermonde);
    // Since d/ddelta = 1 / h d/dxi, the derivative at the evaluation point (xi = 0) is y^(d)(x_e) = d! / h^d c_d, which
    // is linear in the samples, y^(d)(x_e) = sum_i w_i y_i, with the weights w = d! / h^d Q_1 R_1^-T P^T e_d
    // = d! / h^d Q [R_1^-T P^T e_d; 0], computed from right to left in place
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> weights = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>::Zero(offsets.size());
    weights.head(polynomial_degree() + 1) = vandermonde_qr.matrixR()
                                                    .topLeftCorner(polynomial_degree() + 1, polynomial_degree() + 1)
                                                    .template triangularView<Eigen::Upper>()
                                                    .transpose()
                                                    .solve(vandermonde_qr.colsPermutation().transpose() *
                                                            Eigen::Matrix<Scalar, Eigen::Dynamic, 1>::Unit(
                                                                    polynomial_degree() + 1, derivative_order()));
    weights.applyOnTheLeft(vandermonde_qr.householderQ());
    weights *= factorial<Scalar>(derivative_order()) / static_cast<Scalar>(std::pow(scale, derivative_order()));
    return weights;
}

template<std::floating_point Scalar_>
template<typename Samples, IsEuclideanType Sample>
    requires(cppbox::ForwardRangeOf<Samples, Sample>)
inline Sample LocalPolynomialDifferentiator<Scalar_>::weighted_sum(
        const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& weights, const Samples& samples) const {
    check_range_size(samples);
    // Weighted sum, starting from the first term since a dynamic-size zero would need the sample size
    const auto first_sample = std::ranges::begin(samples);
    return std::inner_product(std::next(weights.begin()), weights.end(), std::next(first_sample),
            Sample(weights[0] * (*first_sample)));
}

template<std::floating_point Scalar_>
NonUniformLocalPolynomialDifferentiator<Scalar_>::NonUniformLocalPolynomialDifferentiator(
        const std::size_t window_size_, const unsigned int polynomial_degree_, const unsigned int derivative_order_,
        const std::size_t evaluation_index_)
    : LocalPolynomialDifferentiator<Scalar_>(window_size_, polynomial_degree_, derivative_order_, evaluation_index_) {}

template<std::floating_point Scalar_>
NonUniformLocalPolynomialDifferentiator<Scalar_>::NonUniformLocalPolynomialDifferentiator(
        const std::size_t window_size_, const unsigned int polynomial_degree_, const unsigned int derivative_order_)
    : NonUniformLocalPolynomialDifferentiator(window_size_, polynomial_degree_, derivative_order_,
              LocalPolynomialDifferentiator<Scalar_>::central_index(window_size_)) {}

template<std::floating_point Scalar_>
template<typename SamplePoints, typename Samples, IsIndependentVariableType SamplePoint, IsEuclideanType Sample>
    requires(cppbox::ForwardRangeOf<SamplePoints, SamplePoint> && cppbox::ForwardRangeOf<Samples, Sample>)
inline Sample NonUniformLocalPolynomialDifferentiator<Scalar_>::differentiate(const SamplePoints& sample_points,
        const Samples& samples) const {
    return this->weighted_sum(weights(sample_points), samples);
}

template<std::floating_point Scalar_>
template<typename SamplePoints, typename Samples, IsIndependentVariableType SamplePoint, IsEuclideanType Sample>
    requires(cppbox::ForwardRangeOf<SamplePoints, SamplePoint> && cppbox::ForwardRangeOf<Samples, Sample>)
inline Sample NonUniformLocalPolynomialDifferentiator<Scalar_>::differentiate(const SamplePoints& sample_points,
        const Samples& samples, const SamplePoint& evaluation_point) const {
    return this->weighted_sum(weights(sample_points, evaluation_point), samples);
}

template<std::floating_point Scalar_>
template<typename SamplePoints, IsIndependentVariableType SamplePoint>
    requires(cppbox::ForwardRangeOf<SamplePoints, SamplePoint>)
inline auto NonUniformLocalPolynomialDifferentiator<Scalar_>::weights(const SamplePoints& sample_points) const
        -> Eigen::Matrix<Scalar, Eigen::Dynamic, 1> {
    // Evaluate at the point of the evaluated sample, x_e = x_s, which exists once the number of points is checked
    this->check_range_size(sample_points);
    return weights(sample_points,
            *std::ranges::next(std::ranges::begin(sample_points),
                    static_cast<std::ranges::range_difference_t<SamplePoints>>(this->evaluation_index())));
}

template<std::floating_point Scalar_>
template<typename SamplePoints, IsIndependentVariableType SamplePoint>
    requires(cppbox::ForwardRangeOf<SamplePoints, SamplePoint>)
inline auto NonUniformLocalPolynomialDifferentiator<Scalar_>::weights(const SamplePoints& sample_points,
        const SamplePoint& evaluation_point) const -> Eigen::Matrix<Scalar, Eigen::Dynamic, 1> {
    this->check_range_size(sample_points);
    // Offsets of the sample points from the evaluation point, x_i - x_e, converted to seconds for time points after
    // their subtraction, so that time points far from their clock's epoch keep their precision
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> offsets(this->window_size());
    std::ranges::transform(sample_points, offsets.begin(), [&evaluation_point](const SamplePoint& sample_point) {
        return cppbox::to_sec<Scalar>(sample_point - evaluation_point);
    });
    // Check that the evaluation point is within the window, with sample points at or before it and at or after it,
    // which a NaN evaluation point fails
    throw_if(!((offsets.array() <= static_cast<Scalar>(0)).any() && (offsets.array() >= static_cast<Scalar>(0)).any()),
            "Evaluation point must be within the window.", std::invalid_argument);
    return this->compute_weights(offsets);
}

template<std::floating_point Scalar_>
UniformLocalPolynomialDifferentiator<Scalar_>::UniformLocalPolynomialDifferentiator(const std::size_t window_size_,
        const unsigned int polynomial_degree_, const unsigned int derivative_order_,
        const std::size_t evaluation_index_, const Scalar sample_spacing_)
    : LocalPolynomialDifferentiator<Scalar_>(window_size_, polynomial_degree_, derivative_order_, evaluation_index_),
      sample_spacing_(sample_spacing_) {
    // Check that the sample spacing is positive and finite, where the check is negated so that NaN fails it
    throw_if(!(sample_spacing() > static_cast<Scalar>(0) && std::isfinite(sample_spacing())),
            "Sample spacing must be positive and finite.", std::invalid_argument);
    // Weights for the offsets of the sample points from that of the evaluated sample, x_i - x_s = (i - s) dx, which are
    // the same for every window
    weights_ = this->compute_weights(
            Eigen::Matrix<Scalar, Eigen::Dynamic, 1>::LinSpaced(static_cast<Eigen::Index>(this->window_size()),
                    -static_cast<Scalar>(this->evaluation_index()),
                    static_cast<Scalar>(this->window_size() - 1 - this->evaluation_index())) *
            sample_spacing());
}

template<std::floating_point Scalar_>
UniformLocalPolynomialDifferentiator<Scalar_>::UniformLocalPolynomialDifferentiator(const std::size_t window_size_,
        const unsigned int polynomial_degree_, const unsigned int derivative_order_, const Scalar sample_spacing_)
    : UniformLocalPolynomialDifferentiator(window_size_, polynomial_degree_, derivative_order_,
              LocalPolynomialDifferentiator<Scalar_>::central_index(window_size_), sample_spacing_) {}

template<std::floating_point Scalar_>
template<typename Samples, IsEuclideanType Sample>
    requires(cppbox::ForwardRangeOf<Samples, Sample>)
inline Sample UniformLocalPolynomialDifferentiator<Scalar_>::differentiate(const Samples& samples) const {
    return this->weighted_sum(weights(), samples);
}

template<std::floating_point Scalar_>
inline auto UniformLocalPolynomialDifferentiator<Scalar_>::sample_spacing() const -> Scalar {
    return sample_spacing_;
}

template<std::floating_point Scalar_>
inline auto UniformLocalPolynomialDifferentiator<Scalar_>::weights() const
        -> const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& {
    return weights_;
}

}

#if !MATHBOX_HEADER_ONLY
namespace math {

extern template class LocalPolynomialDifferentiator<double>;
extern template class NonUniformLocalPolynomialDifferentiator<double>;
extern template class UniformLocalPolynomialDifferentiator<double>;

}
#endif

#endif
