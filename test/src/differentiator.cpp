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
#include <random>
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
