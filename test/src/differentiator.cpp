#include "mathbox/differentiator.hpp"

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cppbox/time.hpp>
#include <cstddef>
#include <deque>
#include <limits>
#include <list>
#include <optional>
#include <random>
#include <span>
#include <stdexcept>
#include <vector>

#include "mathbox/combinatorics.hpp"

/**
 * @brief Check that the weights of a uniform differentiator with unit sample spacing equal a known stencil.
 *
 * @param window_size number of samples in the window
 * @param polynomial_degree degree of the fitted polynomial
 * @param derivative_order order of the derivative
 * @param evaluation_index index of the sample at which the derivative is evaluated, counted from the first sample
 * @param expected expected weights
 */
void check_differentiator_weights(const std::size_t window_size, const unsigned int polynomial_degree,
        const unsigned int derivative_order, const std::size_t evaluation_index, const Eigen::VectorXd& expected) {
    const math::UniformLocalPolynomialDifferentiator<double> differentiator(window_size, polynomial_degree,
            derivative_order, evaluation_index, 1.0);
    ASSERT_EQ(differentiator.weights().size(), expected.size());
    EXPECT_TRUE(differentiator.weights().isApprox(expected, 1.0e-12))
            << "weights: " << differentiator.weights().transpose() << "\nexpected: " << expected.transpose();
}

/**
 * @brief Evaluate the derivative of order \f$ m \f$ of a polynomial \f$ P(t) = \sum_k a_k t^k \f$,
 * \f$ P^{(m)}(t) = \sum_{k \geq m} a_k \frac{k!}{(k - m)!} t^{k - m} \f$.
 *
 * @param coefficients coefficients \f$ a_k \f$ of the polynomial, from the constant term
 * @param t time
 * @param order order of the derivative \f$ m \f$
 * @return double
 */
double polynomial_derivative(const std::vector<double>& coefficients, const double t, const unsigned int order) {
    double derivative{0.0};
    for (unsigned int k = order; k < coefficients.size(); ++k) {
        derivative += coefficients[k] * math::factorial<double>(k) / math::factorial<double>(k - order) *
                      std::pow(t, static_cast<double>(k - order));
    }
    return derivative;
}

/**
 * @brief Whether a range of samples can be differentiated at a range of sample points by the non-uniform
 * differentiator.
 *
 * @tparam SamplePoints
 * @tparam Samples
 */
template<class SamplePoints, class Samples>
concept IsNonUniformlyDifferentiable = requires(
        const math::NonUniformLocalPolynomialDifferentiator<double>& differentiator, const SamplePoints& sample_points,
        const Samples& samples) { differentiator.differentiate(sample_points, samples); };

/**
 * @brief Whether a range of samples can be differentiated by the uniform differentiator.
 *
 * @tparam Samples
 */
template<class Samples>
concept IsUniformlyDifferentiable = requires(const math::UniformLocalPolynomialDifferentiator<double>& differentiator,
        const Samples& samples) { differentiator.differentiate(samples); };

/**
 * @brief Whether samples of a type can be buffered and differentiated by the buffered differentiator.
 *
 * @tparam Sample
 */
template<class Sample>
concept IsBufferedDifferentiable = requires { typename math::BufferedLocalPolynomialDifferentiator<double, Sample>; };

/**
 * @brief Get the point of the first sample of the buffered tests, a time far from its clock's epoch.
 *
 * @return std::chrono::sys_time<std::chrono::nanoseconds>
 */
std::chrono::sys_time<std::chrono::nanoseconds> start_time() {
    return std::chrono::sys_time<std::chrono::nanoseconds>{std::chrono::seconds(1'700'000'000)};
}

/**
 * @brief Get the derivatives of an order of three cubic polynomials in the time since the start time.
 *
 * @param time time of the derivatives
 * @param derivative_order order of the derivatives, where order 0 gives the values of the polynomials
 * @return Eigen::Vector3d
 */
Eigen::Vector3d cubic_derivatives(const std::chrono::sys_time<std::chrono::nanoseconds> time,
        const unsigned int derivative_order) {
    const double t = std::chrono::duration<double>(time - start_time()).count();
    return Eigen::Vector3d{polynomial_derivative({0.1, 0.5, -2.0, 3.0}, t, derivative_order),
            polynomial_derivative({-0.2, 1.5, 0.4, -1.0}, t, derivative_order),
            polynomial_derivative({0.3, -0.8, 5.0, 0.5}, t, derivative_order)};
}

/**
 * @brief Get the points of 60 samples 2.5 ms apart from the start time, except for a given spacing between the 30th
 * and the 31st sample.
 *
 * @param spacing spacing of the 30th and the 31st sample
 * @return std::vector<std::chrono::sys_time<std::chrono::nanoseconds>> points of the 60 samples
 */
std::vector<std::chrono::sys_time<std::chrono::nanoseconds>> points_with_middle_spacing(
        const std::chrono::nanoseconds spacing) {
    std::vector<std::chrono::sys_time<std::chrono::nanoseconds>> points;
    for (int i = 0; i < 60; ++i) {
        points.push_back(start_time() + i * std::chrono::microseconds(2500) +
                         (i < 30 ? std::chrono::nanoseconds::zero() : spacing - std::chrono::microseconds(2500)));
    }
    return points;
}

/**
 * @brief Get a buffered differentiator over 11 samples by fitted cubics, with the cubics sampled at some points added.
 *
 * @param points points of the samples, increasing
 * @param max_sample_spacing maximum sample spacing of the differentiator in seconds, if any
 * @param evaluation_index index of the differentiated sample in its window, by default the centre
 * @return math::BufferedLocalPolynomialDifferentiator<std::chrono::sys_time<std::chrono::nanoseconds>, Eigen::Vector3d>
 */
math::BufferedLocalPolynomialDifferentiator<std::chrono::sys_time<std::chrono::nanoseconds>, Eigen::Vector3d>
buffered_cubics(const std::span<const std::chrono::sys_time<std::chrono::nanoseconds>> points,
        const std::optional<double> max_sample_spacing, const std::size_t evaluation_index = 5) {
    math::BufferedLocalPolynomialDifferentiator<std::chrono::sys_time<std::chrono::nanoseconds>, Eigen::Vector3d>
            differentiator(11, 3, evaluation_index, max_sample_spacing);
    for (const std::chrono::sys_time<std::chrono::nanoseconds> point : points) {
        differentiator.add(point, cubic_derivatives(point, 0));
    }
    return differentiator;
}

/**
 * @brief Check that a buffered differentiator over 11 samples, given the cubics sampled at some points, differentiates
 * every sample exactly, at every order up to 3, also over incomplete windows, since they fit cubics exactly too.
 *
 * @param points points of the samples, increasing, in runs without a gap of at least 4 samples
 * @param max_sample_spacing maximum sample spacing of the differentiator in seconds, if any
 * @param evaluation_index index of the differentiated sample in its window, by default the centre
 * @return std::vector<std::chrono::sys_time<std::chrono::nanoseconds>> points of the samples with complete windows
 */
std::vector<std::chrono::sys_time<std::chrono::nanoseconds>> check_cubic_differentiation(
        const std::span<const std::chrono::sys_time<std::chrono::nanoseconds>> points,
        const std::optional<double> max_sample_spacing, const std::size_t evaluation_index = 5) {
    const auto differentiator = buffered_cubics(points, max_sample_spacing, evaluation_index);
    std::vector<std::chrono::sys_time<std::chrono::nanoseconds>> complete_points;
    for (std::size_t i = 0; i < points.size(); ++i) {
        for (unsigned int derivative_order = 0; derivative_order <= 3; ++derivative_order) {
            const std::optional<Eigen::Vector3d> derivatives =
                    differentiator.differentiate(points[i], derivative_order);
            const Eigen::Vector3d expected = cubic_derivatives(points[i], derivative_order);
            EXPECT_TRUE(derivatives && (*derivatives - expected).norm() <= 1.0e-6 * std::max(1.0, expected.norm()))
                    << "derivative order " << derivative_order << " at sample " << i;
        }
        if (differentiator.is_derivative_complete(points[i])) {
            complete_points.push_back(points[i]);
        }
    }
    return complete_points;
}

TEST(NonUniformLocalPolynomialDifferentiator, finite_differences) {
    // Interpolating polynomials (N = p + 1) give the non-uniform finite difference stencils, e.g. the first and second
    // derivatives at the middle of three samples with spacings h_1 and h_2, [-h_2^2, h_2^2 - h_1^2, h_1^2] / (h_1 h_2
    // (h_1 + h_2)) and [2 h_2, -2 (h_1 + h_2), 2 h_1] / (h_1 h_2 (h_1 + h_2)), here with h_1 = 1 and h_2 = 2
    const std::vector<double> times{0.0, 1.0, 3.0};
    EXPECT_TRUE(math::NonUniformLocalPolynomialDifferentiator<double>(3, 2, 1).weights(times).isApprox(
            Eigen::Vector3d{-2.0 / 3.0, 1.0 / 2.0, 1.0 / 6.0}, 1.0e-12));
    EXPECT_TRUE(math::NonUniformLocalPolynomialDifferentiator<double>(3, 2, 2).weights(times).isApprox(
            Eigen::Vector3d{2.0 / 3.0, -1.0, 1.0 / 3.0}, 1.0e-12));
}

TEST(NonUniformLocalPolynomialDifferentiator, uniform_times) {
    // Uniform sample times, as time points or in seconds, give the weights of the uniform differentiator
    const std::chrono::nanoseconds sample_period{2'500'000};
    for (const std::size_t window_size : std::vector<std::size_t>{2, 3, 5, 8, 11}) {
        std::vector<std::chrono::time_point<std::chrono::steady_clock, std::chrono::nanoseconds>> time_points;
        std::vector<double> times;
        for (std::size_t i = 0; i < window_size; ++i) {
            time_points.emplace_back(std::chrono::seconds(100) + static_cast<int>(i) * sample_period);
            times.push_back(0.3 + static_cast<double>(i) * cppbox::to_sec(sample_period));
        }
        for (unsigned int polynomial_degree = 0; polynomial_degree < std::min<std::size_t>(window_size, 5);
                ++polynomial_degree) {
            for (unsigned int derivative_order = 0; derivative_order <= polynomial_degree; ++derivative_order) {
                for (const std::size_t evaluation_index : {std::size_t{0}, (window_size - 1) / 2, window_size - 1}) {
                    const Eigen::VectorXd expected = math::UniformLocalPolynomialDifferentiator<double>(window_size,
                            polynomial_degree, derivative_order, evaluation_index, cppbox::to_sec(sample_period))
                                                             .weights();
                    const math::NonUniformLocalPolynomialDifferentiator<double> differentiator(window_size,
                            polynomial_degree, derivative_order, evaluation_index);
                    EXPECT_TRUE(differentiator.weights(time_points).isApprox(expected, 1.0e-12))
                            << "N = " << window_size << ", p = " << polynomial_degree << ", d = " << derivative_order
                            << ", s = " << evaluation_index;
                    EXPECT_TRUE(differentiator.weights(times).isApprox(expected, 1.0e-12))
                            << "N = " << window_size << ", p = " << polynomial_degree << ", d = " << derivative_order
                            << ", s = " << evaluation_index;
                }
            }
        }
    }
}

TEST(NonUniformLocalPolynomialDifferentiator, centred) {
    const math::NonUniformLocalPolynomialDifferentiator<double> differentiator(7, 3, 2);
    EXPECT_EQ(differentiator.window_size(), 7);
    EXPECT_EQ(differentiator.polynomial_degree(), 3);
    EXPECT_EQ(differentiator.derivative_order(), 2);
    EXPECT_EQ(differentiator.evaluation_index(), 3);
    const std::vector<double> times{0.0, 0.011, 0.019, 0.032, 0.04, 0.048, 0.061};
    EXPECT_TRUE(differentiator.weights(times).isApprox(
            math::NonUniformLocalPolynomialDifferentiator<double>(7, 3, 2, 3).weights(times)));
}

TEST(NonUniformLocalPolynomialDifferentiator, time_points) {
    // Time points far from their clock's epoch give the weights and derivatives of their times since the first sample,
    // since they are subtracted before their conversion to seconds, whereas seconds since the epoch are rounded to
    // about 0.24 microseconds, which changes the weights
    const math::NonUniformLocalPolynomialDifferentiator<double> differentiator(7, 3, 2);
    const std::chrono::time_point<std::chrono::system_clock, std::chrono::nanoseconds> start_time{
            std::chrono::seconds(1'700'000'000)};
    std::vector<std::chrono::time_point<std::chrono::system_clock, std::chrono::nanoseconds>> times;
    std::vector<double> seconds_since_start;
    std::vector<double> samples;
    for (const int nanoseconds_since_start : {0, 2'500'137, 4'999'871, 7'500'412, 10'000'029, 12'499'655, 15'000'318}) {
        const double seconds = 1.0e-9 * nanoseconds_since_start;
        times.push_back(start_time + std::chrono::nanoseconds(nanoseconds_since_start));
        seconds_since_start.push_back(seconds);
        samples.push_back(std::sin(100.0 * seconds));
    }
    EXPECT_TRUE(differentiator.weights(times).isApprox(differentiator.weights(seconds_since_start), 1.0e-12));
    EXPECT_FALSE(differentiator.weights(cppbox::to_secs(times))
                    .isApprox(differentiator.weights(seconds_since_start), 1.0e-6));
    const double expected = differentiator.differentiate(seconds_since_start, samples);
    EXPECT_NEAR(differentiator.differentiate(times, samples), expected, 1.0e-12 * std::abs(expected));
    // Evaluation at a time point between samples matches evaluation at its seconds since the start
    const double expected_between_samples = differentiator.differentiate(seconds_since_start, samples, 0.006);
    EXPECT_NEAR(differentiator.differentiate(times, samples, start_time + std::chrono::nanoseconds(6'000'000)),
            expected_between_samples, 1.0e-12 * std::abs(expected_between_samples));
}

TEST(NonUniformLocalPolynomialDifferentiator, exact_for_polynomials) {
    // Derivatives of every order up to the degree p of random polynomials P(t) = sum_k a_k t^k, sampled at jittered
    // times, at the first, a third, middle and last samples, are exact up to round-off
    std::mt19937 generator(0);
    std::uniform_real_distribution<double> distribution(-1.0, 1.0);
    const double start_time{0.3};
    const double sample_period{0.01};
    for (const std::size_t window_size : std::vector<std::size_t>{2, 3, 5, 8, 11, 21}) {
        // Times jittered by up to 0.4 sample periods, so strictly increasing
        std::vector<double> times;
        for (std::size_t i = 0; i < window_size; ++i) {
            times.push_back(start_time + (static_cast<double>(i) + 0.4 * distribution(generator)) * sample_period);
        }
        for (unsigned int polynomial_degree = 0; polynomial_degree < std::min<std::size_t>(window_size, 5);
                ++polynomial_degree) {
            std::vector<double> coefficients(polynomial_degree + 1);
            for (double& coefficient : coefficients) {
                coefficient = distribution(generator);
            }
            std::vector<double> samples;
            for (const double time : times) {
                samples.push_back(polynomial_derivative(coefficients, time, 0));
            }
            for (unsigned int derivative_order = 0; derivative_order <= polynomial_degree; ++derivative_order) {
                for (const std::size_t evaluation_index :
                        {std::size_t{0}, (window_size - 1) / 3, (window_size - 1) / 2, window_size - 1}) {
                    const math::NonUniformLocalPolynomialDifferentiator<double> differentiator(window_size,
                            polynomial_degree, derivative_order, evaluation_index);
                    const double expected =
                            polynomial_derivative(coefficients, times[evaluation_index], derivative_order);
                    EXPECT_NEAR(differentiator.differentiate(times, samples), expected,
                            1.0e-6 * std::max(1.0, std::abs(expected)))
                            << "N = " << window_size << ", p = " << polynomial_degree << ", d = " << derivative_order
                            << ", s = " << evaluation_index;
                }
            }
        }
    }
}

TEST(NonUniformLocalPolynomialDifferentiator, evaluation_point) {
    // Evaluation at the point of the evaluated sample gives its weights, and evaluation at random points within the
    // window, between the samples, is exact for polynomials of degree up to p
    const std::vector<double> times{0.0, 0.011, 0.019, 0.032, 0.04, 0.048, 0.061};
    for (const std::size_t evaluation_index : {std::size_t{0}, std::size_t{3}, std::size_t{6}}) {
        const math::NonUniformLocalPolynomialDifferentiator<double> differentiator(7, 3, 1, evaluation_index);
        EXPECT_TRUE(differentiator.weights(times, times[evaluation_index]).isApprox(differentiator.weights(times)));
    }
    const std::vector<double> coefficients{0.3, -0.7, 0.5, 0.9};
    std::vector<double> samples;
    for (const double time : times) {
        samples.push_back(polynomial_derivative(coefficients, time, 0));
    }
    std::mt19937 generator(0);
    std::uniform_real_distribution<double> distribution(times.front(), times.back());
    for (int i = 0; i < 20; ++i) {
        const double evaluation_point = distribution(generator);
        for (unsigned int derivative_order = 0; derivative_order <= 3; ++derivative_order) {
            const double expected = polynomial_derivative(coefficients, evaluation_point, derivative_order);
            EXPECT_NEAR(math::NonUniformLocalPolynomialDifferentiator<double>(7, 3, derivative_order)
                                .differentiate(times, samples, evaluation_point),
                    expected, 1.0e-6 * std::max(1.0, std::abs(expected)))
                    << "x_e = " << evaluation_point << ", d = " << derivative_order;
        }
    }
}

TEST(NonUniformLocalPolynomialDifferentiator, vector_samples) {
    // Vector samples are differentiated component-wise
    const math::NonUniformLocalPolynomialDifferentiator<double> differentiator(7, 3, 1);
    const std::vector<double> times{0.0, 0.011, 0.019, 0.032, 0.04, 0.048, 0.061};
    std::vector<Eigen::Vector3d> samples;
    std::vector<double> first_components;
    for (const double t : times) {
        samples.emplace_back(std::sin(t), t * t, 2.0 - t);
        first_components.push_back(std::sin(t));
    }
    const Eigen::Vector3d derivative = differentiator.differentiate(times, samples);
    EXPECT_NEAR(derivative[0], differentiator.differentiate(times, first_components), 1.0e-14);
    EXPECT_NEAR(derivative[0], std::cos(0.032), 1.0e-5);
    EXPECT_NEAR(derivative[1], 0.064, 1.0e-12);
    EXPECT_NEAR(derivative[2], -1.0, 1.0e-12);
    // Dynamic-size samples give the derivative of the fixed-size ones
    const std::vector<Eigen::VectorXd> dynamic_samples(samples.cbegin(), samples.cend());
    EXPECT_TRUE(differentiator.differentiate(times, dynamic_samples).isApprox(derivative, 1.0e-14));
}

TEST(NonUniformLocalPolynomialDifferentiator, ranges) {
    // Times and samples in other forward ranges than vectors, e.g. deques for sliding windows and lists, give the same
    // derivative
    const math::NonUniformLocalPolynomialDifferentiator<double> differentiator(5, 2, 1);
    const std::vector<double> times{0.0, 0.011, 0.019, 0.032, 0.04};
    const std::vector<double> samples{1.0, 1.2, 0.9, 1.4, 1.1};
    const double derivative = differentiator.differentiate(times, samples);
    EXPECT_DOUBLE_EQ(differentiator.differentiate(std::deque<double>(times.cbegin(), times.cend()),
                             std::deque<double>(samples.cbegin(), samples.cend())),
            derivative);
    EXPECT_DOUBLE_EQ(differentiator.differentiate(std::list<double>(times.cbegin(), times.cend()),
                             std::list<double>(samples.cbegin(), samples.cend())),
            derivative);
}

TEST(NonUniformLocalPolynomialDifferentiator, invalid_arguments) {
    // Window size not greater than the polynomial degree, derivative order above the polynomial degree, evaluation
    // index after the window, and an empty or even centred window
    EXPECT_THROW(math::NonUniformLocalPolynomialDifferentiator<double>(3, 3, 1, 1), std::invalid_argument);
    EXPECT_THROW(math::NonUniformLocalPolynomialDifferentiator<double>(5, 2, 3, 2), std::invalid_argument);
    EXPECT_THROW(math::NonUniformLocalPolynomialDifferentiator<double>(5, 2, 1, 5), std::invalid_argument);
    EXPECT_THROW(math::NonUniformLocalPolynomialDifferentiator<double>(0, 0, 0), std::invalid_argument);
    EXPECT_THROW(math::NonUniformLocalPolynomialDifferentiator<double>(6, 2, 1), std::invalid_argument);
    // Numbers of sample points or samples different from the window size, and sample points not strictly increasing or
    // not finite
    const math::NonUniformLocalPolynomialDifferentiator<double> differentiator(3, 2, 1);
    const std::vector<double> samples(3, 0.0);
    EXPECT_THROW(differentiator.differentiate(std::vector<double>{0.0, 1.0}, samples), std::invalid_argument);
    EXPECT_THROW(differentiator.differentiate(std::vector<double>{0.0, 1.0, 2.0}, std::vector<double>(4, 0.0)),
            std::invalid_argument);
    EXPECT_THROW(differentiator.differentiate(std::vector<double>{0.0, 1.0, 1.0}, samples), std::invalid_argument);
    EXPECT_THROW(differentiator.differentiate(std::vector<double>{0.0, 2.0, 1.0}, samples), std::invalid_argument);
    EXPECT_THROW(differentiator.differentiate(std::vector<double>{std::numeric_limits<double>::quiet_NaN(), 0.0, 1.0},
                         samples),
            std::invalid_argument);
    EXPECT_THROW(differentiator.differentiate(std::vector<double>{0.0, 1.0, std::numeric_limits<double>::infinity()},
                         samples),
            std::invalid_argument);
    // Evaluation points outside the window or NaN
    const std::vector<double> times{0.0, 1.0, 2.0};
    EXPECT_THROW(differentiator.differentiate(times, samples, -0.5), std::invalid_argument);
    EXPECT_THROW(differentiator.differentiate(times, samples, 2.5), std::invalid_argument);
    EXPECT_THROW(differentiator.differentiate(times, samples, std::numeric_limits<double>::quiet_NaN()),
            std::invalid_argument);
    // Integer samples are rejected at compile time, since their weighted sum would be truncated, and so are integer
    // sample points, which are neither floating-point values nor time points
    static_assert(IsNonUniformlyDifferentiable<std::vector<double>, std::vector<double>> &&
                  !IsNonUniformlyDifferentiable<std::vector<double>, std::vector<int>> &&
                  !IsNonUniformlyDifferentiable<std::vector<int>, std::vector<double>>);
}

TEST(UniformLocalPolynomialDifferentiator, finite_differences) {
    // Interpolating polynomials (N = p + 1) give the finite difference stencils: the sample itself for a single sample,
    // backward and forward first, central first and second, second-order backward first, and five-point central first
    // differences
    check_differentiator_weights(1, 0, 0, 0, Eigen::VectorXd::Ones(1));
    check_differentiator_weights(2, 1, 1, 1, Eigen::Vector2d{-1.0, 1.0});
    check_differentiator_weights(2, 1, 1, 0, Eigen::Vector2d{-1.0, 1.0});
    check_differentiator_weights(3, 2, 1, 1, Eigen::Vector3d{-0.5, 0.0, 0.5});
    check_differentiator_weights(3, 2, 2, 1, Eigen::Vector3d{1.0, -2.0, 1.0});
    check_differentiator_weights(3, 2, 1, 2, Eigen::Vector3d{0.5, -2.0, 1.5});
    check_differentiator_weights(5, 4, 1, 2, Eigen::Vector<double, 5>{1.0, -8.0, 0.0, 8.0, -1.0} / 12.0);
}

TEST(UniformLocalPolynomialDifferentiator, savitzky_golay) {
    // Least squares polynomials over centred windows give the Savitzky-Golay coefficients, where the odd (even)
    // derivatives of polynomials of degrees 2k - 1 and 2k (2k and 2k + 1) are equal by symmetry
    check_differentiator_weights(5, 2, 0, 2, Eigen::Vector<double, 5>{-3.0, 12.0, 17.0, 12.0, -3.0} / 35.0);
    check_differentiator_weights(5, 2, 1, 2, Eigen::Vector<double, 5>{-2.0, -1.0, 0.0, 1.0, 2.0} / 10.0);
    check_differentiator_weights(5, 3, 1, 2, Eigen::Vector<double, 5>{1.0, -8.0, 0.0, 8.0, -1.0} / 12.0);
    check_differentiator_weights(5, 2, 2, 2, Eigen::Vector<double, 5>{2.0, -1.0, -2.0, -1.0, 2.0} / 7.0);
    check_differentiator_weights(7, 2, 0, 3, Eigen::Vector<double, 7>{-2.0, 3.0, 6.0, 7.0, 6.0, 3.0, -2.0} / 21.0);
    check_differentiator_weights(7, 2, 1, 3, Eigen::Vector<double, 7>{-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0} / 28.0);
    check_differentiator_weights(7, 2, 2, 3, Eigen::Vector<double, 7>{5.0, 0.0, -3.0, -4.0, -3.0, 0.0, 5.0} / 42.0);
    check_differentiator_weights(7, 3, 2, 3, Eigen::Vector<double, 7>{5.0, 0.0, -3.0, -4.0, -3.0, 0.0, 5.0} / 42.0);
}

TEST(UniformLocalPolynomialDifferentiator, centred) {
    const math::UniformLocalPolynomialDifferentiator<double> differentiator(7, 3, 2, 0.01);
    EXPECT_EQ(differentiator.window_size(), 7);
    EXPECT_EQ(differentiator.polynomial_degree(), 3);
    EXPECT_EQ(differentiator.derivative_order(), 2);
    EXPECT_EQ(differentiator.evaluation_index(), 3);
    EXPECT_EQ(differentiator.sample_spacing(), 0.01);
    EXPECT_TRUE(differentiator.weights().isApprox(
            math::UniformLocalPolynomialDifferentiator<double>(7, 3, 2, 3, 0.01).weights()));
}

TEST(UniformLocalPolynomialDifferentiator, sample_spacing) {
    // The weights scale with 1 / dx^d
    for (unsigned int derivative_order = 0; derivative_order <= 3; ++derivative_order) {
        const math::UniformLocalPolynomialDifferentiator<double> unit(9, 3, derivative_order, 2, 1.0);
        const math::UniformLocalPolynomialDifferentiator<double> scaled(9, 3, derivative_order, 2, 0.0025);
        EXPECT_TRUE(scaled.weights().isApprox(unit.weights() / std::pow(0.0025, derivative_order)));
    }
}

TEST(UniformLocalPolynomialDifferentiator, exact_for_polynomials) {
    // Derivatives of every order up to the degree p of random polynomials P(t) = sum_k a_k t^k, at the first, a third,
    // middle and last samples, are exact up to round-off
    std::mt19937 generator(0);
    std::uniform_real_distribution<double> distribution(-1.0, 1.0);
    const double start_time{0.3};
    const double sample_spacing{0.01};
    for (const std::size_t window_size : std::vector<std::size_t>{2, 3, 5, 8, 11, 21}) {
        for (unsigned int polynomial_degree = 0; polynomial_degree < std::min<std::size_t>(window_size, 5);
                ++polynomial_degree) {
            std::vector<double> coefficients(polynomial_degree + 1);
            for (double& coefficient : coefficients) {
                coefficient = distribution(generator);
            }
            std::vector<double> samples;
            for (std::size_t i = 0; i < window_size; ++i) {
                samples.push_back(
                        polynomial_derivative(coefficients, start_time + static_cast<double>(i) * sample_spacing, 0));
            }
            for (unsigned int derivative_order = 0; derivative_order <= polynomial_degree; ++derivative_order) {
                for (const std::size_t evaluation_index :
                        {std::size_t{0}, (window_size - 1) / 3, (window_size - 1) / 2, window_size - 1}) {
                    const math::UniformLocalPolynomialDifferentiator<double> differentiator(window_size,
                            polynomial_degree, derivative_order, evaluation_index, sample_spacing);
                    const double expected = polynomial_derivative(coefficients,
                            start_time + static_cast<double>(evaluation_index) * sample_spacing, derivative_order);
                    EXPECT_NEAR(differentiator.differentiate(samples), expected,
                            1.0e-6 * std::max(1.0, std::abs(expected)))
                            << "N = " << window_size << ", p = " << polynomial_degree << ", d = " << derivative_order
                            << ", s = " << evaluation_index;
                }
            }
        }
    }
}

TEST(UniformLocalPolynomialDifferentiator, vector_samples) {
    // Vector samples are differentiated component-wise
    const math::UniformLocalPolynomialDifferentiator<double> differentiator(7, 3, 1, 0.01);
    std::vector<Eigen::Vector3d> samples;
    std::vector<double> first_components;
    for (int i = 0; i < 7; ++i) {
        const double t = 0.01 * i;
        samples.emplace_back(std::sin(t), t * t, 2.0 - t);
        first_components.push_back(std::sin(t));
    }
    const Eigen::Vector3d derivative = differentiator.differentiate(samples);
    EXPECT_NEAR(derivative[0], differentiator.differentiate(first_components), 1.0e-14);
    EXPECT_NEAR(derivative[0], std::cos(0.03), 1.0e-5);
    EXPECT_NEAR(derivative[1], 0.06, 1.0e-12);
    EXPECT_NEAR(derivative[2], -1.0, 1.0e-12);
}

TEST(UniformLocalPolynomialDifferentiator, invalid_arguments) {
    // Window size not greater than the polynomial degree, derivative order above the polynomial degree, evaluation
    // index after the window, sample spacing not positive, NaN or infinite, and an empty or even centred window
    EXPECT_THROW(math::UniformLocalPolynomialDifferentiator<double>(3, 3, 1, 1, 0.01), std::invalid_argument);
    EXPECT_THROW(math::UniformLocalPolynomialDifferentiator<double>(5, 2, 3, 2, 0.01), std::invalid_argument);
    EXPECT_THROW(math::UniformLocalPolynomialDifferentiator<double>(5, 2, 1, 5, 0.01), std::invalid_argument);
    EXPECT_THROW(math::UniformLocalPolynomialDifferentiator<double>(5, 2, 1, 2, 0.0), std::invalid_argument);
    EXPECT_THROW(
            math::UniformLocalPolynomialDifferentiator<double>(5, 2, 1, 2, std::numeric_limits<double>::quiet_NaN()),
            std::invalid_argument);
    EXPECT_THROW(
            math::UniformLocalPolynomialDifferentiator<double>(5, 2, 1, 2, std::numeric_limits<double>::infinity()),
            std::invalid_argument);
    EXPECT_THROW(math::UniformLocalPolynomialDifferentiator<double>(0, 0, 0, 0.01), std::invalid_argument);
    EXPECT_THROW(math::UniformLocalPolynomialDifferentiator<double>(6, 2, 1, 0.01), std::invalid_argument);
    // Number of samples different from the window size
    const math::UniformLocalPolynomialDifferentiator<double> differentiator(5, 2, 1, 0.01);
    EXPECT_THROW(differentiator.differentiate(std::vector<double>(4, 0.0)), std::invalid_argument);
    // Integer samples are rejected at compile time, since their weighted sum would be truncated
    static_assert(IsUniformlyDifferentiable<std::vector<double>> && !IsUniformlyDifferentiable<std::vector<int>>);
}

TEST(BufferedLocalPolynomialDifferentiator, construction) {
    const math::BufferedLocalPolynomialDifferentiator<std::chrono::sys_time<std::chrono::nanoseconds>, Eigen::Vector3d>
            differentiator(11, 3, 0.0125);
    EXPECT_EQ(differentiator.window_size(), 11);
    EXPECT_EQ(differentiator.polynomial_degree(), 3);
    EXPECT_EQ(differentiator.evaluation_index(), 5);
    EXPECT_EQ(differentiator.max_sample_spacing(), 0.0125);
    EXPECT_EQ(differentiator.size(), 0);
    EXPECT_FALSE(
            (math::BufferedLocalPolynomialDifferentiator<double, double>(5, 2, std::nullopt).max_sample_spacing()));
    // With an evaluation index, e.g. at the end of an even window
    const math::BufferedLocalPolynomialDifferentiator<double, double> last_sample_differentiator(10, 3, 9, 1.0);
    EXPECT_EQ(last_sample_differentiator.window_size(), 10);
    EXPECT_EQ(last_sample_differentiator.evaluation_index(), 9);
    EXPECT_EQ(last_sample_differentiator.max_sample_spacing(), 1.0);
}

TEST(BufferedLocalPolynomialDifferentiator, exact_for_cubics) {
    // Samples of a 400 Hz loop with up to 0.3 ms of jitter, with single missing samples and three missing in a row, all
    // within the maximum spacing of 12.5 ms: all samples are differentiated exactly, and all but the first and last
    // five over complete windows
    std::vector<std::chrono::sys_time<std::chrono::nanoseconds>> points;
    for (int cycle = 0; cycle < 400; ++cycle) {
        if (cycle % 37 == 10 || cycle % 53 == 20 || (cycle >= 200 && cycle < 203)) {
            continue;
        }
        points.push_back(start_time() + std::chrono::round<std::chrono::nanoseconds>(std::chrono::duration<double>(
                                                2.5e-3 * cycle + 0.3e-3 * std::sin(cycle))));
    }
    EXPECT_EQ(check_cubic_differentiation(points, 0.0125),
            std::vector<std::chrono::sys_time<std::chrono::nanoseconds>>(points.begin() + 5, points.end() - 5));

    // Off centre, the samples with s samples before them and 10 - s after them have complete windows
    EXPECT_EQ(check_cubic_differentiation(points, 0.0125, 10),
            std::vector<std::chrono::sys_time<std::chrono::nanoseconds>>(points.begin() + 10, points.end()));
    EXPECT_EQ(check_cubic_differentiation(points, 0.0125, 0),
            std::vector<std::chrono::sys_time<std::chrono::nanoseconds>>(points.begin(), points.end() - 10));
    EXPECT_EQ(check_cubic_differentiation(points, 0.0125, 3),
            std::vector<std::chrono::sys_time<std::chrono::nanoseconds>>(points.begin() + 3, points.end() - 7));
}

TEST(BufferedLocalPolynomialDifferentiator, final_windows) {
    // As samples are added, the window of each sample is final and complete once the fifth sample after it is added,
    // while the windows of the first five samples, shifted to start at the first sample, are final once the eleventh
    // sample is added, and never complete. Every sample has derivatives once four samples are buffered.
    const std::vector<std::chrono::sys_time<std::chrono::nanoseconds>> points =
            points_with_middle_spacing(std::chrono::microseconds(2500));
    auto differentiator = buffered_cubics({}, 0.0125);
    for (std::size_t i = 0; i < 30; ++i) {
        differentiator.add(points[i], cubic_derivatives(points[i], 0));
        for (std::size_t j = 0; j <= i; ++j) {
            EXPECT_EQ(differentiator.is_derivative_final(points[j]), j < 5 ? i >= 10 : j + 5 <= i)
                    << "sample " << j << " of " << i;
            EXPECT_EQ(differentiator.is_derivative_complete(points[j]), j >= 5 && j + 5 <= i)
                    << "sample " << j << " of " << i;
            EXPECT_EQ(differentiator.differentiate(points[j], 1).has_value(), i >= 3) << "sample " << j << " of " << i;
            EXPECT_EQ(differentiator.is_derivative_computable(points[j]), i >= 3) << "sample " << j << " of " << i;
        }
    }
}

TEST(BufferedLocalPolynomialDifferentiator, final_windows_off_centre) {
    // With windows of three samples before and seven after their sample, the window of each sample is final and
    // complete once the seventh sample after it is added, while the windows of the first three samples are final once
    // the eleventh sample is added, and never complete
    const std::vector<std::chrono::sys_time<std::chrono::nanoseconds>> points =
            points_with_middle_spacing(std::chrono::microseconds(2500));
    auto differentiator = buffered_cubics({}, 0.0125, 3);
    for (std::size_t i = 0; i < 30; ++i) {
        differentiator.add(points[i], cubic_derivatives(points[i], 0));
        for (std::size_t j = 0; j <= i; ++j) {
            EXPECT_EQ(differentiator.is_derivative_final(points[j]), j < 3 ? i >= 10 : j + 7 <= i)
                    << "sample " << j << " of " << i;
            EXPECT_EQ(differentiator.is_derivative_complete(points[j]), j >= 3 && j + 7 <= i)
                    << "sample " << j << " of " << i;
        }
    }
}

TEST(BufferedLocalPolynomialDifferentiator, final_windows_ending_at_the_sample) {
    // With windows ending at their sample, the window of each sample from the eleventh on is final and complete as soon
    // as it is added, while the windows of the first ten are final once the eleventh sample is added, and never
    // complete
    const std::vector<std::chrono::sys_time<std::chrono::nanoseconds>> points =
            points_with_middle_spacing(std::chrono::microseconds(2500));
    auto differentiator = buffered_cubics({}, 0.0125, 10);
    for (std::size_t i = 0; i < 30; ++i) {
        differentiator.add(points[i], cubic_derivatives(points[i], 0));
        EXPECT_EQ(differentiator.is_derivative_final_up_to(points[i]), i >= 10) << "sample " << i;
        for (std::size_t j = 0; j <= i; ++j) {
            EXPECT_EQ(differentiator.is_derivative_final(points[j]), i >= 10) << "sample " << j << " of " << i;
            EXPECT_EQ(differentiator.is_derivative_complete(points[j]), j >= 10) << "sample " << j << " of " << i;
        }
    }
}

TEST(BufferedLocalPolynomialDifferentiator, final_windows_up_to) {
    // As samples with jitter, missing samples and a gap of 13 ms after the 40th are added, the windows of all samples
    // up to a point, at a sample or between samples, are final exactly when the window of each of them is, so windows
    // become final in the order of their samples. Without samples up to the point, they are final at once.
    std::vector<std::chrono::sys_time<std::chrono::nanoseconds>> points{start_time()};
    for (int i = 1; i < 80; ++i) {
        points.push_back(
                points.back() + (i == 40 ? std::chrono::microseconds(13000)
                                         : std::chrono::round<std::chrono::nanoseconds>(std::chrono::duration<double>(
                                                   (i % 17 == 0 ? 5.0e-3 : 2.5e-3) + 0.3e-3 * std::sin(i)))));
    }
    auto differentiator = buffered_cubics({}, 0.0125);
    EXPECT_TRUE(differentiator.is_derivative_final_up_to(points[0]));
    for (std::size_t i = 0; i < points.size(); ++i) {
        differentiator.add(points[i], cubic_derivatives(points[i], 0));
        EXPECT_TRUE(differentiator.is_derivative_final_up_to(points[0] - std::chrono::nanoseconds(1)));
        bool expected{true};
        for (std::size_t j = 0; j <= i; ++j) {
            expected = expected && differentiator.is_derivative_final(points[j]);
            EXPECT_EQ(differentiator.is_derivative_final_up_to(points[j]), expected) << "sample " << j << " of " << i;
            EXPECT_EQ(differentiator.is_derivative_final_up_to(points[j] + std::chrono::nanoseconds(1)), expected)
                    << "after sample " << j << " of " << i;
        }
    }
}

TEST(BufferedLocalPolynomialDifferentiator, gaps) {
    // A spacing of 13 ms after the 30th sample, larger than the maximum of 12.5 ms, so a gap: all samples are
    // differentiated exactly, and those with five samples on each side and no gap between them over complete windows
    const std::vector<std::chrono::sys_time<std::chrono::nanoseconds>> points =
            points_with_middle_spacing(std::chrono::microseconds(13000));
    std::vector<std::chrono::sys_time<std::chrono::nanoseconds>> expected_points(points.begin() + 5,
            points.begin() + 25);
    expected_points.insert(expected_points.end(), points.begin() + 35, points.begin() + 55);
    EXPECT_EQ(check_cubic_differentiation(points, 0.0125), expected_points);
    // With windows ending at their sample, those with ten samples before them and no gap between them
    std::vector<std::chrono::sys_time<std::chrono::nanoseconds>> expected_last_points(points.begin() + 10,
            points.begin() + 30);
    expected_last_points.insert(expected_last_points.end(), points.begin() + 40, points.end());
    EXPECT_EQ(check_cubic_differentiation(points, 0.0125, 10), expected_last_points);

    // The windows of the last five samples before the gap are final once the first sample after it is added, shifted to
    // end before the gap, so incomplete
    auto differentiator = buffered_cubics(std::span(points).first(30), 0.0125);
    for (std::size_t i = 25; i < 30; ++i) {
        EXPECT_FALSE(differentiator.is_derivative_final(points[i]));
    }
    differentiator.add(points[30], cubic_derivatives(points[30], 0));
    for (std::size_t i = 25; i < 30; ++i) {
        EXPECT_TRUE(differentiator.is_derivative_final(points[i]));
        EXPECT_FALSE(differentiator.is_derivative_complete(points[i]));
    }

    // After the gap, the windows of the first six samples, shifted to start after the gap, are final once eleven
    // samples after the gap are buffered, and only the sixth is complete
    for (std::size_t i = 31; i < 40; ++i) {
        differentiator.add(points[i], cubic_derivatives(points[i], 0));
    }
    EXPECT_FALSE(differentiator.is_derivative_final(points[30]));
    differentiator.add(points[40], cubic_derivatives(points[40], 0));
    for (std::size_t i = 30; i < 36; ++i) {
        EXPECT_TRUE(differentiator.is_derivative_final(points[i]));
        EXPECT_EQ(differentiator.is_derivative_complete(points[i]), i == 35);
    }
    EXPECT_FALSE(differentiator.is_derivative_final(points[36]));
}

TEST(BufferedLocalPolynomialDifferentiator, spacing_of_the_maximum) {
    // A spacing of exactly the maximum of 12.5 ms is spanned, so all but the first and last five samples have complete
    // windows, while a spacing of 1 ns more is a gap
    const std::vector<std::chrono::sys_time<std::chrono::nanoseconds>> spanned =
            points_with_middle_spacing(std::chrono::microseconds(12500));
    EXPECT_EQ(check_cubic_differentiation(spanned, 0.0125),
            std::vector<std::chrono::sys_time<std::chrono::nanoseconds>>(spanned.begin() + 5, spanned.end() - 5));
    const std::vector<std::chrono::sys_time<std::chrono::nanoseconds>> split =
            points_with_middle_spacing(std::chrono::microseconds(12500) + std::chrono::nanoseconds(1));
    EXPECT_EQ(check_cubic_differentiation(split, 0.0125).size(), 40);

    // So the window of a sample that spans the spacing of the maximum is not final until half a window of samples is
    // buffered after it
    const auto differentiator = buffered_cubics(std::span(spanned).first(32), 0.0125);
    EXPECT_TRUE(differentiator.is_derivative_final(spanned[26]));
    EXPECT_FALSE(differentiator.is_derivative_final(spanned[27]));
}

TEST(BufferedLocalPolynomialDifferentiator, no_maximum_sample_spacing) {
    // Without a maximum sample spacing, a spacing of 1 s is spanned, so all but the first and last five samples have
    // complete windows
    const std::vector<std::chrono::sys_time<std::chrono::nanoseconds>> points =
            points_with_middle_spacing(std::chrono::seconds(1));
    EXPECT_EQ(check_cubic_differentiation(points, std::nullopt),
            std::vector<std::chrono::sys_time<std::chrono::nanoseconds>>(points.begin() + 5, points.end() - 5));

    // So the windows of the samples before the spacing are not final until half a window of samples is buffered after
    // them
    const auto differentiator = buffered_cubics(std::span(points).first(31), std::nullopt);
    EXPECT_TRUE(differentiator.is_derivative_final(points[25]));
    EXPECT_FALSE(differentiator.is_derivative_final(points[26]));
}

TEST(BufferedLocalPolynomialDifferentiator, samples_as_they_arrive) {
    // Samples with jitter and missing samples, and with spacings of the maximum, just above it and well above it, are
    // added one by one, and each sample is differentiated once its window is final, after which the samples before it
    // are removed: the derivatives equal those of a buffer of all the samples, and only the newest samples remain
    // pending
    std::vector<std::chrono::sys_time<std::chrono::nanoseconds>> points{start_time()};
    for (int i = 1; i < 200; ++i) {
        const std::chrono::nanoseconds spacing = [i]() -> std::chrono::nanoseconds {
            if (i == 60) {
                return std::chrono::microseconds(12500);
            }
            if (i == 110) {
                return std::chrono::microseconds(12500) + std::chrono::nanoseconds(1);
            }
            if (i == 150) {
                return std::chrono::milliseconds(50);
            }
            return std::chrono::round<std::chrono::nanoseconds>(
                    std::chrono::duration<double>((i % 23 == 0 ? 5.0e-3 : 2.5e-3) + 0.3e-3 * std::sin(i)));
        }();
        points.push_back(points.back() + spacing);
    }
    const auto full_buffer = buffered_cubics(points, 0.0125);
    auto differentiator = buffered_cubics({}, 0.0125);
    std::size_t next = 0;
    std::size_t complete = 0;
    for (std::size_t i = 0; i < points.size(); ++i) {
        differentiator.add(points[i], cubic_derivatives(points[i], 0));
        for (; next <= i && differentiator.is_derivative_final(points[next]); ++next) {
            for (unsigned int derivative_order = 1; derivative_order <= 2; ++derivative_order) {
                EXPECT_EQ(differentiator.differentiate(points[next], derivative_order),
                        full_buffer.differentiate(points[next], derivative_order))
                        << "sample " << next << ", derivative order " << derivative_order;
            }
            EXPECT_EQ(differentiator.is_derivative_complete(points[next]),
                    full_buffer.is_derivative_complete(points[next]))
                    << "sample " << next;
            complete += differentiator.is_derivative_complete(points[next]);
            differentiator.remove_before(points[next]);
        }
    }
    // Gaps after the 110th and the 150th sample, so all windows are complete except the five on each side of each gap
    // and at each end, and the last five remain pending
    EXPECT_EQ(complete, 170);
    EXPECT_EQ(next, 195);
    for (std::size_t i = next; i < points.size(); ++i) {
        EXPECT_FALSE(differentiator.is_derivative_final(points[i]));
        EXPECT_FALSE(full_buffer.is_derivative_complete(points[i]));
    }
}

TEST(BufferedLocalPolynomialDifferentiator, remove_before) {
    // Removing the samples before a sample keeps the last ten of them, the window size less one, which a window shifted
    // to fit its run may need
    const std::vector<std::chrono::sys_time<std::chrono::nanoseconds>> points =
            points_with_middle_spacing(std::chrono::microseconds(2500));
    auto differentiator = buffered_cubics(std::span(points).first(30), 0.0125);
    const std::optional<Eigen::Vector3d> second_derivatives = differentiator.differentiate(points[15], 2);
    differentiator.remove_before(points[15]);
    EXPECT_EQ(differentiator.size(), 25);
    EXPECT_TRUE(differentiator.is_derivative_complete(points[15]));
    EXPECT_EQ(differentiator.differentiate(points[15], 2), second_derivatives);
    // Before a point between samples, it keeps the last ten samples before the next sample
    differentiator.remove_before(points[20] - std::chrono::nanoseconds(1));
    EXPECT_EQ(differentiator.size(), 20);
    // Before a point after the samples, it keeps the last ten of them
    differentiator.remove_before(points[29] + std::chrono::seconds(1));
    EXPECT_EQ(differentiator.size(), 10);
    // Before a point within ten samples of the first buffered sample, it keeps them all
    differentiator.remove_before(points[27]);
    EXPECT_EQ(differentiator.size(), 10);

    // With windows ending at their sample, the window of the sample is unchanged, while those of the kept samples
    // before it are incomplete, as if near the first buffered sample
    auto last_sample_differentiator = buffered_cubics(std::span(points).first(30), 0.0125, 10);
    const std::optional<Eigen::Vector3d> last_sample_derivatives =
            last_sample_differentiator.differentiate(points[15], 1);
    last_sample_differentiator.remove_before(points[15]);
    EXPECT_EQ(last_sample_differentiator.size(), 25);
    EXPECT_EQ(last_sample_differentiator.differentiate(points[15], 1), last_sample_derivatives);
    EXPECT_FALSE(last_sample_differentiator.is_derivative_complete(points[14]));
}

TEST(BufferedLocalPolynomialDifferentiator, newest_samples) {
    // The windows of the newest five samples are neither final nor complete, but shifted to end at the newest sample,
    // they give derivatives, exact for cubics, e.g. at the end of the data
    const std::vector<std::chrono::sys_time<std::chrono::nanoseconds>> points =
            points_with_middle_spacing(std::chrono::microseconds(2500));
    const auto differentiator = buffered_cubics(std::span(points).first(30), 0.0125);
    EXPECT_TRUE(differentiator.is_derivative_final_up_to(points[24]));
    EXPECT_FALSE(differentiator.is_derivative_final_up_to(points[25]));
    for (std::size_t i = 25; i < 30; ++i) {
        EXPECT_FALSE(differentiator.is_derivative_final(points[i])) << "sample " << i;
        EXPECT_FALSE(differentiator.is_derivative_complete(points[i])) << "sample " << i;
        const std::optional<Eigen::Vector3d> second_derivatives = differentiator.differentiate(points[i], 2);
        EXPECT_TRUE(second_derivatives && second_derivatives->isApprox(cubic_derivatives(points[i], 2), 1.0e-6))
                << "sample " << i;
    }
}

TEST(BufferedLocalPolynomialDifferentiator, short_runs) {
    // Runs of 20, 6, 3 and 15 samples between gaps of 20 ms: the run of 6 gives derivatives over a window shortened to
    // it, exact for cubics and final once the gap after it is buffered, while the run of 3 has fewer samples than a
    // cubic needs, so no derivatives
    std::vector<std::chrono::sys_time<std::chrono::nanoseconds>> points;
    std::chrono::sys_time<std::chrono::nanoseconds> point = start_time();
    for (const int run : {20, 6, 3, 15}) {
        for (int i = 0; i < run; ++i) {
            points.push_back(point);
            point += std::chrono::microseconds(2500);
        }
        point += std::chrono::milliseconds(20);
    }
    auto differentiator = buffered_cubics(std::span(points).first(26), 0.0125);
    EXPECT_FALSE(differentiator.is_derivative_final(points[22]));
    differentiator.add(points[26], cubic_derivatives(points[26], 0));
    for (std::size_t i = 20; i < 26; ++i) {
        EXPECT_TRUE(differentiator.is_derivative_final(points[i])) << "sample " << i;
        EXPECT_FALSE(differentiator.is_derivative_complete(points[i])) << "sample " << i;
        EXPECT_TRUE(differentiator.is_derivative_computable(points[i])) << "sample " << i;
        const std::optional<Eigen::Vector3d> derivatives = differentiator.differentiate(points[i], 1);
        EXPECT_TRUE(derivatives && derivatives->isApprox(cubic_derivatives(points[i], 1), 1.0e-6)) << "sample " << i;
    }
    for (std::size_t i = 27; i < points.size(); ++i) {
        differentiator.add(points[i], cubic_derivatives(points[i], 0));
    }
    for (std::size_t i = 26; i < 29; ++i) {
        EXPECT_TRUE(differentiator.is_derivative_final(points[i])) << "sample " << i;
        EXPECT_FALSE(differentiator.is_derivative_computable(points[i])) << "sample " << i;
        EXPECT_FALSE(differentiator.differentiate(points[i], 1)) << "sample " << i;
    }
}

TEST(BufferedLocalPolynomialDifferentiator, scalar_samples) {
    // A quadratic sampled at jittered floating-point points, differentiated by quadratics fitted over 5 samples
    math::BufferedLocalPolynomialDifferentiator<double, double> differentiator(5, 2, 0.1);
    const std::vector<double> points{0.0, 0.051, 0.098, 0.152, 0.2, 0.249, 0.301};
    for (const double point : points) {
        differentiator.add(point, 1.0 + 2.0 * point - 3.0 * point * point);
    }
    for (std::size_t i = 0; i < points.size(); ++i) {
        const std::optional<double> derivative = differentiator.differentiate(points[i], 1);
        const std::optional<double> second_derivative = differentiator.differentiate(points[i], 2);
        EXPECT_EQ(differentiator.is_derivative_complete(points[i]), i >= 2 && i < 5);
        ASSERT_TRUE(derivative && second_derivative);
        EXPECT_NEAR(*derivative, 2.0 - 6.0 * points[i], 1.0e-9);
        EXPECT_NEAR(*second_derivative, -6.0, 1.0e-9);
    }
}

TEST(BufferedLocalPolynomialDifferentiator, matrix_samples) {
    // Dynamic-size matrices of quadratics sampled at jittered points, differentiated entry by entry by quadratics
    // fitted over 5 samples, while samples of another number of rows or columns are rejected
    math::BufferedLocalPolynomialDifferentiator<double, Eigen::MatrixXd> differentiator(5, 2, 0.1);
    const std::vector<double> points{0.0, 0.051, 0.098, 0.152, 0.2};
    for (const double t : points) {
        Eigen::MatrixXd sample(2, 2);
        sample << 1.0 + t, 2.0 - t * t, 3.0 * t * t, -t + 0.5 * t * t;
        differentiator.add(t, sample);
    }
    const double t = points[2];
    Eigen::MatrixXd expected_derivatives(2, 2);
    expected_derivatives << 1.0, -2.0 * t, 6.0 * t, -1.0 + t;
    Eigen::MatrixXd expected_second_derivatives(2, 2);
    expected_second_derivatives << 0.0, -2.0, 6.0, 1.0;
    const std::optional<Eigen::MatrixXd> derivatives = differentiator.differentiate(t, 1);
    const std::optional<Eigen::MatrixXd> second_derivatives = differentiator.differentiate(t, 2);
    ASSERT_TRUE(derivatives && second_derivatives);
    EXPECT_TRUE(derivatives->isApprox(expected_derivatives, 1.0e-9)) << *derivatives;
    EXPECT_TRUE(second_derivatives->isApprox(expected_second_derivatives, 1.0e-9)) << *second_derivatives;
    EXPECT_THROW(differentiator.add(0.25, Eigen::MatrixXd::Zero(2, 3)), std::invalid_argument);
    EXPECT_THROW(differentiator.add(0.25, Eigen::MatrixXd::Zero(3, 2)), std::invalid_argument);
    EXPECT_EQ(differentiator.size(), 5);
}

TEST(BufferedLocalPolynomialDifferentiator, invalid_arguments) {
    // An even centred window, an evaluation index outside the window, a window not greater than the polynomial degree,
    // and a maximum sample spacing that is not positive, infinite or NaN
    EXPECT_THROW((math::BufferedLocalPolynomialDifferentiator<double, double>(10, 3, 1.0)), std::invalid_argument);
    EXPECT_THROW((math::BufferedLocalPolynomialDifferentiator<double, double>(11, 3, 11, 1.0)), std::invalid_argument);
    EXPECT_THROW((math::BufferedLocalPolynomialDifferentiator<double, double>(3, 3, 1.0)), std::invalid_argument);
    EXPECT_THROW((math::BufferedLocalPolynomialDifferentiator<double, double>(5, 2, 0.0)), std::invalid_argument);
    EXPECT_THROW((math::BufferedLocalPolynomialDifferentiator<double, double>(5, 2, -1.0)), std::invalid_argument);
    EXPECT_THROW((math::BufferedLocalPolynomialDifferentiator<double, double>(5, 2,
                         std::numeric_limits<double>::infinity())),
            std::invalid_argument);
    EXPECT_THROW((math::BufferedLocalPolynomialDifferentiator<double, double>(5, 2,
                         std::numeric_limits<double>::quiet_NaN())),
            std::invalid_argument);
    // A NaN or infinite sample point is rejected, even as the first, which would otherwise block the buffer
    math::BufferedLocalPolynomialDifferentiator<double, double> scalar_differentiator(5, 2, 1.0);
    EXPECT_THROW(scalar_differentiator.add(std::numeric_limits<double>::quiet_NaN(), 0.0), std::invalid_argument);
    EXPECT_THROW(scalar_differentiator.add(std::numeric_limits<double>::infinity(), 0.0), std::invalid_argument);
    EXPECT_EQ(scalar_differentiator.size(), 0);
    // A sample at or before the last sample, or of another size, is rejected without changing the buffer
    math::BufferedLocalPolynomialDifferentiator<double, Eigen::VectorXd> differentiator(5, 2, 1.0);
    differentiator.add(0.0, Eigen::VectorXd::Zero(3));
    EXPECT_THROW(differentiator.add(0.0, Eigen::VectorXd::Zero(3)), std::invalid_argument);
    EXPECT_THROW(differentiator.add(-1.0, Eigen::VectorXd::Zero(3)), std::invalid_argument);
    EXPECT_THROW(differentiator.add(1.0, Eigen::VectorXd::Zero(4)), std::invalid_argument);
    EXPECT_EQ(differentiator.size(), 1);
    // A derivative order above the polynomial degree, and points without a buffered sample
    EXPECT_THROW(differentiator.differentiate(0.0, 3), std::invalid_argument);
    EXPECT_THROW(differentiator.differentiate(0.5, 1), std::invalid_argument);
    EXPECT_THROW(differentiator.differentiate(-0.5, 1), std::invalid_argument);
    EXPECT_THROW(differentiator.is_derivative_final(0.5), std::invalid_argument);
    EXPECT_THROW(differentiator.is_derivative_final(-0.5), std::invalid_argument);
    EXPECT_THROW(differentiator.is_derivative_complete(0.5), std::invalid_argument);
    EXPECT_THROW(differentiator.is_derivative_computable(0.5), std::invalid_argument);
    // Samples of integer scalars are rejected at compile time, since their weighted sum would be truncated
    static_assert(IsBufferedDifferentiable<double> && IsBufferedDifferentiable<Eigen::Vector3d> &&
                  !IsBufferedDifferentiable<int> && !IsBufferedDifferentiable<Eigen::Vector3i>);
}
