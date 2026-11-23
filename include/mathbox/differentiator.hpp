#ifndef MATHBOX_DIFFERENTIATOR_HPP
#define MATHBOX_DIFFERENTIATOR_HPP

#include <Eigen/Core>
#include <concepts>
#include <cppbox/ranges.hpp>
#include <cstddef>
#include <deque>
#include <optional>
#include <ranges>
#include <vector>

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
 * - Evaluating at an end of the window gives a noisier derivative than inside it, unless \f$ d = p \f$. The last
 *   sample needs no later samples, while the central sample needs half a window of them, which delays its derivative
 *   by half a window when the samples arrive one by one.
 * - Samples within half a window of the start or end of the data, or of a gap in it, have no centred window that
 *   stays within the data. They can be evaluated off centre, e.g. the first sample after a gap with \f$ s = 0 \f$, at
 *   the cost of noisier derivatives near the ends of the window, as above, while
 *   BufferedLocalPolynomialDifferentiator gives them no derivative.
 *
 * This base holds and checks \f$ N \f$, \f$ p \f$, \f$ d \f$ and \f$ s \f$. The derived classes provide the weights:
 * UniformLocalPolynomialDifferentiator computes them once for uniformly spaced samples, and
 * NonUniformLocalPolynomialDifferentiator computes them for each window from its sample points.
 * BufferedLocalPolynomialDifferentiator applies the latter to samples that arrive one by one.
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
     * @brief Check that a range has one element for each sample of the window.
     *
     * @tparam Range forward range
     * @param range range of the samples or of their points
     */
    template<std::ranges::forward_range Range>
    void check_range_size(const Range& range) const;

    /**
     * @brief Compute the weights of the samples from the offsets of their points from the evaluation point,
     * \f$ \delta_i = x_i - x_e \f$, e.g. the point \f$ x_s \f$ of the evaluated sample.
     *
     * The weights are found in three steps:
     * 1. The offsets are scaled to \f$ \xi_i = \delta_i / h \in [-1, 1] \f$ with \f$ h = \max_i |\delta_i| \f$, or
     *    \f$ h = 1 \f$ for a single sample at the evaluation point, so that their powers are in \f$ [-1, 1] \f$ for any
     *    unit of \f$ x \f$.
     * 2. The polynomial \f$ y(\xi) = \sum_{j=0}^{p} c_j \xi^j \f$ is fitted to the samples \f$ y \f$ with the
     *    column-pivoted QR decomposition \f$ A P = Q R \f$ of the Vandermonde matrix \f$ A_{ij} = \xi_i^j \f$, whose
     *    leading \f$ (p + 1) \times (p + 1) \f$ block \f$ R_1 \f$ of \f$ R \f$ is invertible since the offsets are
     *    distinct, which gives the coefficients \f$ c = P R_1^{-1} Q_1^T y \f$.
     * 3. Since \f$ d/d\delta = h^{-1} d/d\xi \f$, the derivative is \f$ y^{(d)}(x_e) = d! \, c_d / h^d \f$, which is
     *    the weighted sum of the samples with the weights
     *    \f[
     *        w = \frac{d!}{h^d} Q \begin{bmatrix} R_1^{-T} P^T e_d \\ 0 \end{bmatrix},
     *    \f]
     *    where \f$ e_d \f$ is the unit vector of the coefficient \f$ c_d \f$.
     *
     * The round-off error of the weights grows with the polynomial degree, and it is largest for an evaluation point at
     * an end of the window.
     *
     * @param offsets offsets \f$ \delta_i = x_i - x_e \f$ of the points of the \f$ N \f$ samples from the evaluation
     * point, finite and strictly increasing
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
 * The derivative is evaluated at the evaluated sample or, given an evaluation point \f$ x_e \f$, at any point within
 * the window, e.g. the time of another measurement.
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
 * // Or at a time within the window, e.g. that of another measurement
 * const std::chrono::steady_clock::time_point time = stamps[3] + std::chrono::milliseconds(1);
 * const Eigen::Matrix<double, 12, 1> derivatives_at_time = differentiator.differentiate(stamps, values, time);
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
     * @brief Differentiate a window of samples at their points, at an evaluation point within the window,
     * \f$ y^{(d)}(x_e) = \sum_{i=0}^{N-1} w_i y_i \f$.
     *
     * @tparam SamplePoints forward range of sample points
     * @tparam Samples forward range of samples
     * @tparam SamplePoint type of the sample points and the evaluation point, a floating-point scalar or a
     * `std::chrono::time_point`
     * @tparam Sample type of the samples, deduced from them, a floating-point scalar or an Eigen matrix, differentiated
     * component-wise
     * @param sample_points points \f$ x_i \f$ of the \f$ N \f$ samples, strictly increasing
     * @param samples window of \f$ N \f$ samples, in the order of their points
     * @param evaluation_point point \f$ x_e \f$ at which the derivative is evaluated, within the window,
     * \f$ x_0 \leq x_e \leq x_{N-1} \f$
     * @return Sample derivative at the evaluation point
     */
    template<typename SamplePoints, typename Samples,
            IsIndependentVariableType SamplePoint = std::ranges::range_value_t<SamplePoints>,
            IsEuclideanType Sample = std::ranges::range_value_t<Samples>>
        requires(cppbox::ForwardRangeOf<SamplePoints, SamplePoint> && cppbox::ForwardRangeOf<Samples, Sample>)
    Sample differentiate(const SamplePoints& sample_points, const Samples& samples,
            const SamplePoint& evaluation_point) const;

    /**
     * @brief Compute the weights \f$ w_i \f$ of the samples for their points, at the evaluated sample, in the order of
     * the samples.
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

    /**
     * @brief Compute the weights \f$ w_i \f$ of the samples for their points, at an evaluation point within the window,
     * in the order of the samples.
     *
     * @tparam SamplePoints forward range of sample points
     * @tparam SamplePoint type of the sample points and the evaluation point, a floating-point scalar or a
     * `std::chrono::time_point`
     * @param sample_points points \f$ x_i \f$ of the \f$ N \f$ samples, strictly increasing
     * @param evaluation_point point \f$ x_e \f$ at which the derivative is evaluated, within the window,
     * \f$ x_0 \leq x_e \leq x_{N-1} \f$
     * @return Eigen::Matrix<Scalar, Eigen::Dynamic, 1>
     */
    template<typename SamplePoints, IsIndependentVariableType SamplePoint = std::ranges::range_value_t<SamplePoints>>
        requires(cppbox::ForwardRangeOf<SamplePoints, SamplePoint>)
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> weights(const SamplePoints& sample_points,
            const SamplePoint& evaluation_point) const;
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
     * @param sample_spacing_ spacing \f$ \Delta x \f$ of the sample points, positive and finite
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
     * @param sample_spacing_ spacing \f$ \Delta x \f$ of the sample points, positive and finite
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

/**
 * @brief Differentiator of samples that arrive one by one at increasing points, which buffers them and differentiates a
 * buffered sample by a local polynomial fit over its window of samples (see NonUniformLocalPolynomialDifferentiator).
 *
 * The samples are added as they arrive, e.g. measurements as they are received:
 * 1. Window: a sample is differentiated over its window, the \f$ N \f$ samples from \f$ s \f$ before it to
 *    \f$ N - 1 - s \f$ after it, with the evaluation index \f$ s \f$. The centred window, \f$ s = (N - 1) / 2 \f$,
 *    delays the derivatives by half a window, while the window ending at the sample, \f$ s = N - 1 \f$, needs no
 *    later samples, so it has no delay but noisier derivatives (see LocalPolynomialDifferentiator).
 *    Unlike the latency, the derivatives have no lag, i.e. no delay of the signal as a causal low-pass filter has:
 *    the fit is exact for polynomials of degree up to \f$ p \f$, so their derivatives are at the sample itself for any
 *    \f$ s \f$, and the centred window of uniformly spaced samples has symmetric weights, so no phase shift at any
 *    frequency, only attenuation of the fast changes. Off centre, the weights are asymmetric, so content faster than
 *    the window fits well is shifted in phase, and amplified rather than attenuated near an end of the window.
 * 2. Edges and gaps: with a maximum sample spacing, a window must not span a gap, i.e. two consecutive samples
 *    spaced by more than it. Where the window of a sample does not fit in its run of buffered samples without a gap,
 *    i.e. near the first buffered sample, next to a gap, or near the newest sample until later samples arrive, the
 *    window is shifted to fit the run, and shortened to the run if the run has fewer than \f$ N \f$ samples. The
 *    derivatives from such a window are incomplete (see `is_derivative_complete()`), so less accurate, as point 1
 *    describes off centre, and not computable below \f$ p + 1 \f$ samples (see `is_derivative_computable()`).
 * 3. Final derivatives: the derivatives of a sample are final once adding samples can no longer change them, or its
 *    window: once it is buffered as the window of \f$ N \f$ samples that the run allows, or once the run is ended
 *    by a gap (see `is_derivative_final()`). A sample processed as samples arrive should be differentiated once its
 *    derivatives are final, and they used, or not, depending on whether they are complete. Derivatives become final
 *    in the order of their samples, so those of all samples up to a point are final once those of the last one are
 *    (see `is_derivative_final_up_to()`).
 * 4. Removal: the buffer keeps every added sample until removed, so samples that no later window needs can be removed
 *    explicitly (see `remove_before()`), e.g. after differentiating each sample in turn, to bound its memory and the
 *    cost of finding a sample in it.
 *
 * Choosing the maximum sample spacing:
 * - A window that spans missing samples is fitted at the actual sample points, so its derivatives remain valid, but
 *   they become less accurate as the spacing grows, the higher derivatives most. The maximum sample spacing is the
 *   largest spacing that a window may span, and without one, windows span any spacing.
 * - Like the sample spacing of UniformLocalPolynomialDifferentiator, it is a floating-point value, in seconds for time
 *   points, whose spacings are converted to seconds after their subtraction.
 * - To accept up to \f$ k \f$ missing samples in a row, it lies between \f$ k + 1 \f$ and \f$ k + 2 \f$ sample periods,
 *   e.g. \f$ k + 1.5 \f$ periods so that jitter does not turn \f$ k \f$ missing samples into a gap.
 * - Samples near a gap get an incomplete window shifted to avoid the gap, whose derivatives other than of order
 *   \f$ d = p \f$ are much less accurate than from their own window across a short gap (see point 1), so a short
 *   run of missing samples is better spanned than treated as a gap.
 *
 * Example: the first and second derivatives of 12 values sampled at about 400 Hz, from fitted cubics over 11 samples,
 * with up to 4 missing samples in a row accepted by a maximum spacing of 5.5 sample periods,
 * @code{.cpp}
 * math::BufferedLocalPolynomialDifferentiator<std::chrono::steady_clock::time_point, Eigen::Matrix<double, 12, 1>>
 *         differentiator(11, 3, 0.01375);
 * // For each sample, as it arrives
 * differentiator.add(stamp, values);
 * // For a buffered sample, e.g. that of the next measurement to process, once its derivatives are final
 * if (differentiator.is_derivative_final(sample_stamp)) {
 *     const std::optional<Eigen::Matrix<double, 12, 1>> derivatives = differentiator.differentiate(sample_stamp, 1);
 *     // Only from a complete window, since second derivatives from an incomplete one are much less accurate
 *     std::optional<Eigen::Matrix<double, 12, 1>> second_derivatives;
 *     if (differentiator.is_derivative_complete(sample_stamp)) {
 *         second_derivatives = differentiator.differentiate(sample_stamp, 2);
 *     }
 *     differentiator.remove_before(sample_stamp);
 * }
 * @endcode
 * where the derivatives of a sample are final once the samples up to 5 samples (12.5 ms) after it are buffered. With
 * `differentiator(11, 3, 10, 0.01375)` instead, each sample is differentiated over the 11 samples up to it, so its
 * derivatives are final as soon as it is added.
 *
 * @tparam SamplePoint_ type of the sample points, a floating-point scalar or a `std::chrono::time_point`
 * @tparam Sample_ type of the samples, a floating-point scalar or an Eigen matrix of floating-point scalars,
 * differentiated component-wise
 */
template<IsIndependentVariableType SamplePoint_, IsEuclideanType Sample_>
    requires(std::floating_point<typename MathTypeTraits<Sample_>::Scalar>)
class BufferedLocalPolynomialDifferentiator {
public:
    /**
     * @brief Type of the sample points.
     *
     */
    using SamplePoint = SamplePoint_;

    /**
     * @brief Type of the samples.
     *
     */
    using Sample = Sample_;

    /**
     * @brief Floating-point scalar type of the samples.
     *
     */
    using Scalar = typename MathTypeTraits<Sample>::Scalar;

    /**
     * @brief Construct an empty buffer.
     *
     * @param window_size_ number of samples in the window \f$ N \f$, greater than the polynomial degree
     * @param polynomial_degree_ degree of the fitted polynomials \f$ p \f$, the highest order of the derivatives
     * @param evaluation_index_ index \f$ s \f$ of the differentiated sample in its window, counted from the first
     * sample of the window, in \f$ \{0, \dots, N - 1\} \f$ (e.g. \f$ N - 1 \f$ for windows ending at the sample)
     * @param max_sample_spacing_ largest spacing of consecutive sample points that a window may span, positive and
     * finite, in seconds for time points, or none for windows that span any spacing
     */
    explicit BufferedLocalPolynomialDifferentiator(const std::size_t window_size_,
            const unsigned int polynomial_degree_, const std::size_t evaluation_index_,
            const std::optional<Scalar> max_sample_spacing_);

    /**
     * @brief Construct an empty buffer that differentiates each sample over the window centred on it,
     * \f$ s = (N - 1) / 2 \f$, where the window size must be odd so that its centre is a sample.
     *
     * @param window_size_ number of samples in the window \f$ N \f$, odd and greater than the polynomial degree
     * @param polynomial_degree_ degree of the fitted polynomials \f$ p \f$, the highest order of the derivatives
     * @param max_sample_spacing_ largest spacing of consecutive sample points that a window may span, positive and
     * finite, in seconds for time points, or none for windows that span any spacing
     */
    explicit BufferedLocalPolynomialDifferentiator(const std::size_t window_size_,
            const unsigned int polynomial_degree_, const std::optional<Scalar> max_sample_spacing_);

    /**
     * @brief Add a sample after the buffered samples.
     *
     * @param sample_point point of the sample, after the points of the buffered samples
     * @param sample sample, of the size of the buffered samples
     */
    void add(const SamplePoint& sample_point, const Sample& sample);

    /**
     * @brief Differentiate the buffered sample at a point over its window, \f$ y^{(d)}(x_s) = \sum_i w_i y_i \f$,
     * where an incomplete window (see `is_derivative_complete()`) is shifted or shortened to fit the run of samples
     * without a gap that contains the sample.
     *
     * @param sample_point point of a buffered sample
     * @param derivative_order order of the derivative \f$ d \f$, at most the polynomial degree
     * @return std::optional<Sample> derivative at the sample, per second to the power \f$ d \f$ for time points, or
     * `std::nullopt` if they are not computable (see `is_derivative_computable()`)
     */
    std::optional<Sample> differentiate(const SamplePoint& sample_point, const unsigned int derivative_order) const;

    /**
     * @brief Get the index \f$ s \f$ of the differentiated sample in its window, counted from the first sample of the
     * window, i.e. the number of samples before it in its window.
     *
     * @return std::size_t
     */
    std::size_t evaluation_index() const;

    /**
     * @brief Check if the derivatives of the buffered sample at a point are complete, i.e. if its window is the window
     * of \f$ N \f$ samples from \f$ s \f$ before the sample to \f$ N - 1 - s \f$ after it, all buffered and spanning
     * no gap, so that its derivatives are as accurate as configured. Complete derivatives are final.
     *
     * @param sample_point point of a buffered sample
     * @return true
     * @return false
     */
    bool is_derivative_complete(const SamplePoint& sample_point) const;

    /**
     * @brief Check if the derivatives of the buffered sample at a point can be computed (see `differentiate()`),
     * i.e. if its window, shifted or shortened to fit its run of samples without a gap, has more than \f$ p \f$
     * samples.
     *
     * @param sample_point point of a buffered sample
     * @return true
     * @return false
     */
    bool is_derivative_computable(const SamplePoint& sample_point) const;

    /**
     * @brief Check if the derivatives of the buffered sample at a point are final, i.e. if they (see
     * `differentiate()`), and its window, no longer change as samples are added: once the window of \f$ N \f$
     * samples that its run of samples without a gap allows is buffered, or once the run is ended by a gap. They may be
     * incomplete (see `is_derivative_complete()`) or not computable (see `is_derivative_computable()`).
     *
     * @param sample_point point of a buffered sample
     * @return true
     * @return false
     */
    bool is_derivative_final(const SamplePoint& sample_point) const;

    /**
     * @brief Check if the derivatives of all the buffered samples at or before a point are final (see
     * `is_derivative_final()`), e.g. to wait until all samples up to a point can be differentiated.
     *
     * @param sample_point point, of a buffered sample or not
     * @return true if they are, including if no sample at or before the point is buffered
     * @return false
     */
    bool is_derivative_final_up_to(const SamplePoint& sample_point) const;

    /**
     * @brief Get the largest spacing of consecutive sample points that a window may span, if any, in seconds for time
     * points.
     *
     * @return std::optional<Scalar>
     */
    std::optional<Scalar> max_sample_spacing() const;

    /**
     * @brief Get the degree of the fitted polynomials \f$ p \f$.
     *
     * @return unsigned int
     */
    unsigned int polynomial_degree() const;

    /**
     * @brief Remove the buffered samples before a point that the windows of the samples at or after it do not need,
     * i.e. all but the last \f$ N - 1 \f$ of them, which a window shifted to fit its run may need.
     *
     * The samples before the point should then no longer be queried: the removed ones throw, and the windows of the
     * kept ones lack samples, so they are differentiated as if near the first buffered sample.
     *
     * @param sample_point point, e.g. that of the last differentiated sample when the samples are differentiated in
     * increasing order
     */
    void remove_before(const SamplePoint& sample_point);

    /**
     * @brief Get the number of buffered samples.
     *
     * @return std::size_t
     */
    std::size_t size() const;

    /**
     * @brief Get the number of samples in the window \f$ N \f$.
     *
     * @return std::size_t
     */
    std::size_t window_size() const;

private:
    /**
     * @brief Find the index of the buffered sample at exactly a point, from the oldest buffered sample.
     *
     * @param sample_point point
     * @return std::optional<std::size_t> index of the sample, or `std::nullopt` if none is at the point
     */
    std::optional<std::size_t> find_index(const SamplePoint& sample_point) const;

    /**
     * @brief Window over which the derivatives of a buffered sample are computed: its configured window, or, if
     * incomplete, that window shifted or shortened to fit the sample's run of samples without a gap.
     *
     */
    struct DerivativeWindow {
        // Index of the first sample of the window
        std::size_t first;
        // Index after the last sample of the window
        std::size_t end;
        // Whether the window, and so the derivatives, no longer change as samples are added
        bool final;
    };

    /**
     * @brief Get the derivative window of the buffered sample at an index (see DerivativeWindow): its configured window
     * if it is buffered and spans no gap, and otherwise the window of at most \f$ N \f$ samples nearest to it within
     * its run of samples without a gap (see point 2 of the class description).
     *
     * @param index index of the sample
     * @return DerivativeWindow
     */
    DerivativeWindow derivative_window(const std::size_t index) const;

    /**
     * @brief Get the number of samples after the differentiated sample in its window, \f$ N - 1 - s \f$, e.g.
     * \f$ (N - 1) / 2 \f$, half a window, for the centred window.
     *
     * @return std::size_t
     */
    std::size_t num_samples_after() const;

    /**
     * @brief Check if the buffered samples with indices from `first` to `last - 1` span a gap, i.e. if there is a
     * maximum sample spacing and two consecutive ones are spaced by more than it.
     *
     * @param first index of the first sample
     * @param last index after the last sample
     * @return true
     * @return false
     */
    bool spans_gap(const std::size_t first, const std::size_t last) const;

    // Differentiators at the evaluation index, indexed by their derivative order, from 0 to the polynomial degree
    std::vector<NonUniformLocalPolynomialDifferentiator<Scalar>> differentiators_;
    // Largest spacing of consecutive sample points that a window may span, if any
    std::optional<Scalar> max_sample_spacing_;
    // Buffered sample points, increasing
    std::deque<SamplePoint> sample_points_;
    // Buffered samples, in the order of their points
    std::deque<Sample> samples_;
};

}

#include "mathbox/impl/differentiator.hpp"

#endif
