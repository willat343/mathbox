#ifndef MATHBOX_DIFFERENTIATOR_HPP
#define MATHBOX_DIFFERENTIATOR_HPP

#include <Eigen/Core>
#include <concepts>
#include <cppbox/ranges.hpp>
#include <cstddef>
#include <ranges>

#include "mathbox/traits.hpp"

namespace math {

/**
 * @brief Base of the differentiators of sampled values by a local polynomial fit.
 *
 * A window of \f$ N \f$ samples \f$ y_0, \dots, y_{N-1} \f$ is taken at strictly increasing points
 * \f$ x_0, \dots, x_{N-1} \f$ of an independent variable \f$ x \f$, e.g. time. Their derivative of order \f$ d \f$ is
 * estimated in three steps:
 * 1. A polynomial of degree \f$ p \f$ is fitted to the samples in the least squares sense.
 * 2. Its derivative of order \f$ d \f$ is evaluated at the point \f$ x_s \f$ of the sample with index
 *    \f$ s \in \{0, \dots, N - 1\} \f$, the evaluation index.
 * 3. Since the fit is linear in the samples, this derivative is a weighted sum of them,
 *    \f$ y^{(d)}(x_s) = \sum_{i=0}^{N-1} w_i y_i \f$, whose weights only depend on the offsets \f$ x_i - x_s \f$ (see
 *    compute_weights).
 *
 * Choosing the parameters:
 * - The derivative is exact for polynomials of degree up to \f$ p \f$. Finite differences are the case
 *   \f$ N = p + 1 \f$, where the polynomial interpolates the samples.
 * - Windows of more than \f$ p + 1 \f$ samples smooth the derivative, attenuating noise but also fast changes.
 * - Evaluating at a sample at an end of the window gives a noisier derivative than at one inside it, unless
 *   \f$ d = p \f$. The last sample needs no later samples, while the central sample needs half a window of them,
 *   which delays its derivative by half a window when the samples arrive one by one.
 * - Samples within half a window of the start or end of the data, or of a gap in it, have no centred window that
 *   stays within the data, so they are evaluated off centre, e.g. the first sample after a gap with \f$ s = 0 \f$.
 *
 * This base holds and checks \f$ N \f$, \f$ p \f$, \f$ d \f$ and \f$ s \f$. The derived classes provide the weights:
 * UniformLocalPolynomialDifferentiator computes them once for uniformly spaced samples, and
 * NonUniformLocalPolynomialDifferentiator computes them for each window from its sample points.
 *
 * @tparam Scalar_ floating-point scalar type
 */
template<std::floating_point Scalar_>
class LocalPolynomialDifferentiator {
public:
    /**
     * @brief Floating-point scalar type.
     *
     */
    using Scalar = Scalar_;

    /**
     * @brief Get the order of the derivative \f$ d \f$.
     *
     * @return unsigned int
     */
    unsigned int derivative_order() const;

    /**
     * @brief Get the index \f$ s \f$ of the sample at which the derivative is evaluated, counted from the first
     * sample of the window.
     *
     * @return std::size_t
     */
    std::size_t evaluation_index() const;

    /**
     * @brief Get the degree of the fitted polynomial \f$ p \f$.
     *
     * @return unsigned int
     */
    unsigned int polynomial_degree() const;

    /**
     * @brief Get the number of samples in the window \f$ N \f$.
     *
     * @return std::size_t
     */
    std::size_t window_size() const;

protected:
    /**
     * @brief Construct the base, checking the parameters.
     *
     * @param window_size_ number of samples in the window \f$ N \f$, greater than the polynomial degree
     * @param polynomial_degree_ degree of the fitted polynomial \f$ p \f$, at least the derivative order
     * @param derivative_order_ order of the derivative \f$ d \f$, where order 0 smooths the samples
     * @param evaluation_index_ index \f$ s \f$ of the sample at which the derivative is evaluated, counted from the
     * first sample of the window, in \f$ \{0, \dots, N - 1\} \f$
     */
    explicit LocalPolynomialDifferentiator(const std::size_t window_size_, const unsigned int polynomial_degree_,
            const unsigned int derivative_order_, const std::size_t evaluation_index_);

    /**
     * @brief Get the index of the central sample of a window, \f$ s = (N - 1) / 2 \f$, for the centred constructors,
     * checking that the window size is odd so that its centre is a sample.
     *
     * @param window_size number of samples in the window \f$ N \f$, odd
     * @return std::size_t
     */
    static std::size_t central_index(const std::size_t window_size);

    /**
     * @brief Compute the weights of the samples from the offsets of their points from that of the evaluated sample,
     * \f$ \delta_i = x_i - x_s \f$.
     *
     * The weights are found in three steps:
     * 1. The offsets are scaled to \f$ \xi_i = \delta_i / h \in [-1, 1] \f$ with \f$ h = \max_i |\delta_i| \f$, or
     *    \f$ h = 1 \f$ for a single sample, so that their powers are in \f$ [-1, 1] \f$ for any unit of \f$ x \f$.
     * 2. The polynomial \f$ y(\xi) = \sum_{j=0}^{p} c_j \xi^j \f$ is fitted to the samples \f$ y \f$ with the
     *    column-pivoted QR decomposition \f$ A P = Q R \f$ of the Vandermonde matrix \f$ A_{ij} = \xi_i^j \f$, whose
     *    leading \f$ (p + 1) \times (p + 1) \f$ block \f$ R_1 \f$ of \f$ R \f$ is invertible since the offsets are
     *    distinct, which gives the coefficients \f$ c = P R_1^{-1} Q_1^T y \f$.
     * 3. Since \f$ d/d\delta = h^{-1} d/d\xi \f$, the derivative is \f$ y^{(d)}(x_s) = d! \, c_d / h^d \f$, which is
     *    the weighted sum of the samples with the weights
     *    \f[
     *        w = \frac{d!}{h^d} Q \begin{bmatrix} R_1^{-T} P^T e_d \\ 0 \end{bmatrix},
     *    \f]
     *    where \f$ e_d \f$ is the unit vector of the coefficient \f$ c_d \f$.
     *
     * The round-off error of the weights grows with the polynomial degree, and it is largest for an evaluated sample at
     * an end of the window.
     *
     * @param offsets offsets \f$ \delta_i = x_i - x_s \f$ of the points of the \f$ N \f$ samples from that of the
     * evaluated sample, finite and strictly increasing
     * @return Eigen::Matrix<Scalar, Eigen::Dynamic, 1> weights \f$ w_i \f$ of the samples, in the order of the samples
     */
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> compute_weights(
            const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& offsets) const;

    /**
     * @brief Compute the weighted sum of a window of samples, \f$ \sum_{i=0}^{N-1} w_i y_i \f$.
     *
     * @tparam Samples forward range of samples
     * @tparam Sample type of the samples, deduced from them, a floating-point scalar or an Eigen matrix
     * @param weights weights \f$ w_i \f$ of the samples, in the order of the samples
     * @param samples window of \f$ N \f$ samples, in the order of their points
     * @return Sample
     */
    template<typename Samples, IsEuclideanType Sample = std::ranges::range_value_t<Samples>>
        requires(cppbox::ForwardRangeOf<Samples, Sample>)
    Sample weighted_sum(const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& weights, const Samples& samples) const;

private:
    // Number of samples in the window
    std::size_t window_size_;
    // Degree of the fitted polynomial
    unsigned int polynomial_degree_;
    // Order of the derivative
    unsigned int derivative_order_;
    // Index of the sample at which the derivative is evaluated, counted from the first sample of the window
    std::size_t evaluation_index_;
};

/**
 * @brief Differentiator of non-uniformly sampled values by a local polynomial fit (see LocalPolynomialDifferentiator).
 *
 * The weights depend on the sample points of each window, so they are computed for each window by a QR decomposition
 * (see LocalPolynomialDifferentiator::compute_weights), which costs much more than the weighted sum that follows.
 *
 * The sample points are floating-point values or time points:
 * - Time points are subtracted before their conversion to seconds, exactly for integer durations such as nanoseconds,
 *   so time points far from their clock's epoch keep their precision, and the derivative is per second to the power
 *   \f$ d \f$.
 * - Floating-point values carry the round-off of their magnitude, e.g. about 0.24 microseconds for seconds since 1970,
 *   so they should be relative to a nearby origin.
 *
 * Noise in the sample points becomes noise in the weights, and so in the derivative. Samples taken at a uniform rate
 * but stamped with noisy times are therefore differentiated more accurately by UniformLocalPolynomialDifferentiator
 * with the nominal sample spacing.
 *
 * Example: the first derivatives of 12 values at the central one of 7 time stamps, each from a fitted quadratic,
 * @code{.cpp}
 * const math::NonUniformLocalPolynomialDifferentiator<double> differentiator(7, 2, 1);
 * // Time stamps and values of the latest 7 samples, from the oldest
 * std::deque<std::chrono::steady_clock::time_point> stamps;
 * std::deque<Eigen::Matrix<double, 12, 1>> values;
 * const Eigen::Matrix<double, 12, 1> derivatives = differentiator.differentiate(stamps, values);
 * @endcode
 *
 * @tparam Scalar_ floating-point scalar type
 */
template<std::floating_point Scalar_>
class NonUniformLocalPolynomialDifferentiator : public LocalPolynomialDifferentiator<Scalar_> {
public:
    /**
     * @brief Floating-point scalar type.
     *
     */
    using Scalar = Scalar_;

    /**
     * @brief Construct the differentiator.
     *
     * @param window_size_ number of samples in the window \f$ N \f$, greater than the polynomial degree
     * @param polynomial_degree_ degree of the fitted polynomial \f$ p \f$, at least the derivative order
     * @param derivative_order_ order of the derivative \f$ d \f$, where order 0 smooths the samples
     * @param evaluation_index_ index \f$ s \f$ of the sample at which the derivative is evaluated, counted from the
     * first sample of the window, in \f$ \{0, \dots, N - 1\} \f$ (e.g. \f$ (N - 1) / 2 \f$ for the centre of an odd
     * window, \f$ N - 1 \f$ for the last sample)
     */
    explicit NonUniformLocalPolynomialDifferentiator(const std::size_t window_size_,
            const unsigned int polynomial_degree_, const unsigned int derivative_order_,
            const std::size_t evaluation_index_);

    /**
     * @brief Construct the differentiator evaluating at the centre of the window, \f$ s = (N - 1) / 2 \f$, where the
     * window size must be odd so that its centre is a sample.
     *
     * @param window_size_ number of samples in the window \f$ N \f$, odd and greater than the polynomial degree
     * @param polynomial_degree_ degree of the fitted polynomial \f$ p \f$, at least the derivative order
     * @param derivative_order_ order of the derivative \f$ d \f$, where order 0 smooths the samples
     */
    explicit NonUniformLocalPolynomialDifferentiator(const std::size_t window_size_,
            const unsigned int polynomial_degree_, const unsigned int derivative_order_);

    /**
     * @brief Differentiate a window of samples at their points, \f$ y^{(d)}(x_s) = \sum_{i=0}^{N-1} w_i y_i \f$.
     *
     * @tparam SamplePoints forward range of sample points
     * @tparam Samples forward range of samples
     * @tparam SamplePoint type of the sample points, deduced from them, a floating-point scalar or a
     * `std::chrono::time_point`
     * @tparam Sample type of the samples, deduced from them, a floating-point scalar or an Eigen matrix, differentiated
     * component-wise
     * @param sample_points points \f$ x_i \f$ of the \f$ N \f$ samples, strictly increasing
     * @param samples window of \f$ N \f$ samples, in the order of their points
     * @return Sample derivative at the evaluated sample
     */
    template<typename SamplePoints, typename Samples,
            IsIndependentVariableType SamplePoint = std::ranges::range_value_t<SamplePoints>,
            IsEuclideanType Sample = std::ranges::range_value_t<Samples>>
        requires(cppbox::ForwardRangeOf<SamplePoints, SamplePoint> && cppbox::ForwardRangeOf<Samples, Sample>)
    Sample differentiate(const SamplePoints& sample_points, const Samples& samples) const;

    /**
     * @brief Compute the weights \f$ w_i \f$ of the samples for their points, in the order of the samples.
     *
     * @tparam SamplePoints forward range of sample points
     * @tparam SamplePoint type of the sample points, deduced from them, a floating-point scalar or a
     * `std::chrono::time_point`
     * @param sample_points points \f$ x_i \f$ of the \f$ N \f$ samples, strictly increasing
     * @return Eigen::Matrix<Scalar, Eigen::Dynamic, 1>
     */
    template<typename SamplePoints, IsIndependentVariableType SamplePoint = std::ranges::range_value_t<SamplePoints>>
        requires(cppbox::ForwardRangeOf<SamplePoints, SamplePoint>)
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> weights(const SamplePoints& sample_points) const;
};

/**
 * @brief Differentiator of uniformly sampled values by a local polynomial fit (see LocalPolynomialDifferentiator), also
 * known as a Savitzky-Golay filter.
 *
 * The samples are at the points \f$ x_i = x_0 + i \Delta x \f$ with the sample spacing \f$ \Delta x \f$, e.g. the
 * sample period of a signal in time. The weights therefore only depend on \f$ N \f$, \f$ p \f$, \f$ d \f$,
 * \f$ s \f$ and \f$ \Delta x \f$, so they are computed once, on construction.
 *
 * Finite differences are the special case \f$ N = p + 1 \f$, where the polynomial interpolates the samples, with the
 * weights (times \f$ \Delta x^d \f$):
 * - \f$ N = 2, p = 1, d = 1 \f$: backward difference \f$ [-1, 1] \f$ with \f$ s = 1 \f$, forward difference with
 *   \f$ s = 0 \f$,
 * - \f$ N = 3, p = 2, d = 1, s = 1 \f$: central difference \f$ [-\frac{1}{2}, 0, \frac{1}{2}] \f$,
 * - \f$ N = 3, p = 2, d = 2, s = 1 \f$: central second difference \f$ [1, -2, 1] \f$,
 * - \f$ N = 3, p = 2, d = 1, s = 2 \f$: second-order backward difference \f$ [\frac{1}{2}, -2, \frac{3}{2}] \f$,
 * - \f$ N = 5, p = 4, d = 1, s = 2 \f$: five-point central difference \f$ [1, -8, 0, 8, -1] / 12 \f$.
 *
 * Example: the second derivatives of 12 values sampled at 400 Hz, at the central one of 11 samples, each from a fitted
 * cubic,
 * @code{.cpp}
 * const math::UniformLocalPolynomialDifferentiator<double> differentiator(11, 3, 2, 0.0025);
 * // Values of the latest 11 samples, from the oldest
 * std::deque<Eigen::Matrix<double, 12, 1>> values;
 * const Eigen::Matrix<double, 12, 1> second_derivatives = differentiator.differentiate(values);
 * @endcode
 * where the second derivatives are at the sample 5 samples (12.5 ms) before the latest one.
 *
 * @tparam Scalar_ floating-point scalar type
 */
template<std::floating_point Scalar_>
class UniformLocalPolynomialDifferentiator : public LocalPolynomialDifferentiator<Scalar_> {
public:
    /**
     * @brief Floating-point scalar type.
     *
     */
    using Scalar = Scalar_;

    /**
     * @brief Construct the differentiator, computing its weights.
     *
     * @param window_size_ number of samples in the window \f$ N \f$, greater than the polynomial degree
     * @param polynomial_degree_ degree of the fitted polynomial \f$ p \f$, at least the derivative order
     * @param derivative_order_ order of the derivative \f$ d \f$, where order 0 smooths the samples
     * @param evaluation_index_ index \f$ s \f$ of the sample at which the derivative is evaluated, counted from the
     * first sample of the window, in \f$ \{0, \dots, N - 1\} \f$ (e.g. \f$ (N - 1) / 2 \f$ for the centre of an odd
     * window, \f$ N - 1 \f$ for the last sample)
     * @param sample_spacing_ spacing \f$ \Delta x \f$ of the sample points, positive
     */
    explicit UniformLocalPolynomialDifferentiator(const std::size_t window_size_, const unsigned int polynomial_degree_,
            const unsigned int derivative_order_, const std::size_t evaluation_index_, const Scalar sample_spacing_);

    /**
     * @brief Construct the differentiator evaluating at the centre of the window, \f$ s = (N - 1) / 2 \f$, computing
     * its weights, where the window size must be odd so that its centre is a sample.
     *
     * @param window_size_ number of samples in the window \f$ N \f$, odd and greater than the polynomial degree
     * @param polynomial_degree_ degree of the fitted polynomial \f$ p \f$, at least the derivative order
     * @param derivative_order_ order of the derivative \f$ d \f$, where order 0 smooths the samples
     * @param sample_spacing_ spacing \f$ \Delta x \f$ of the sample points, positive
     */
    explicit UniformLocalPolynomialDifferentiator(const std::size_t window_size_, const unsigned int polynomial_degree_,
            const unsigned int derivative_order_, const Scalar sample_spacing_);

    /**
     * @brief Differentiate a window of samples, \f$ y^{(d)}(x_s) = \sum_{i=0}^{N-1} w_i y_i \f$.
     *
     * @tparam Samples forward range of samples
     * @tparam Sample type of the samples, deduced from them, a floating-point scalar or an Eigen matrix, differentiated
     * component-wise
     * @param samples window of \f$ N \f$ samples, in the order of their points
     * @return Sample derivative at the evaluated sample
     */
    template<typename Samples, IsEuclideanType Sample = std::ranges::range_value_t<Samples>>
        requires(cppbox::ForwardRangeOf<Samples, Sample>)
    Sample differentiate(const Samples& samples) const;

    /**
     * @brief Get the spacing \f$ \Delta x \f$ of the sample points.
     *
     * @return Scalar
     */
    Scalar sample_spacing() const;

    /**
     * @brief Get the weights \f$ w_i \f$ of the samples, in the order of the samples.
     *
     * @return const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>&
     */
    const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& weights() const;

private:
    // Spacing of the sample points
    Scalar sample_spacing_;
    // Weights of the samples, in the order of the samples
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> weights_;
};

}

#include "mathbox/impl/differentiator.hpp"

#endif
