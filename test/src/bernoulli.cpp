#include "mathbox/bernoulli.hpp"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <stdexcept>

TEST(bernoulli, function) {
    EXPECT_DOUBLE_EQ(math::bernoulli::function(0.0), 1.0);
    // Arguments in the series and in the closed form of the quadratic coefficient
    for (const double x : {-5.0, -1.5, -0.3, 1.0e-6, 0.3, 0.999, 1.001, 4.0, 20.0}) {
        EXPECT_NEAR(math::bernoulli::function(x), x / std::expm1(x), 1.0e-12 * (1.0 + std::abs(x))) << x;
    }
}

TEST(bernoulli, number) {
    // Known values at compile time, the odd numbers after B_1 being zero
    static_assert(math::bernoulli::number<double>(0) == 1.0 && math::bernoulli::number<double>(1) == -0.5 &&
                  math::bernoulli::number<double>(2) == 1.0 / 6.0 && math::bernoulli::number<double>(3) == 0.0 &&
                  math::bernoulli::number<double>(4) == -1.0 / 30.0 && math::bernoulli::number<double>(29) == 0.0);
    EXPECT_EQ(math::bernoulli::number<double>(12), -691.0 / 2730.0);
    EXPECT_EQ(math::bernoulli::number<double>(20), -174611.0 / 330.0);
    EXPECT_EQ(math::bernoulli::number<double>(30), 8615841276005.0 / 14322.0);
    EXPECT_THROW(math::bernoulli::number<double>(31), std::out_of_range);
}

TEST(bernoulli, quadratic_coefficient_closed_forms) {
    // Real arguments in the series, either side of the change to the closed form, and in the closed form, with the
    // Langevin function L(x) = coth(x) - 1 / x
    for (const double x : {1.0e-3, 0.2, 0.999, 1.001, 2.0, 30.0}) {
        const double expected = (1.0 / std::tanh(0.5 * x) - 2.0 / x) / (2.0 * x);
        EXPECT_NEAR(math::bernoulli::quadratic_coefficient<0>(x * x)[0], expected, 1.0e-9 * expected) << x;
    }
    // Angles (imaginary arguments), likewise, and at pi
    for (const double angle : {1.0e-3, 0.2, 0.999, 1.001, 2.0, M_PI, 5.0}) {
        const double expected = (1.0 - 0.5 * angle / std::tan(0.5 * angle)) / (angle * angle);
        EXPECT_NEAR(math::bernoulli::quadratic_coefficient<0>(-angle * angle)[0], expected, 1.0e-9 * expected) << angle;
    }
    EXPECT_DOUBLE_EQ(math::bernoulli::quadratic_coefficient<0>(0.0)[0], 1.0 / 12.0);
}

TEST(bernoulli, quadratic_coefficient_partial_fractions) {
    // q(z) = sum_n 2 / (4 pi^2 n^2 + z), whose tail after N terms is about 1 / (2 pi^2 N)
    for (const double squared_argument : {-9.0, -0.5, 0.5, 9.0}) {
        double sum = 0.0;
        const int num_terms = 2000000;
        for (int n = num_terms; n >= 1; --n) {
            sum += 2.0 / (4.0 * M_PI * M_PI * n * n + squared_argument);
        }
        EXPECT_NEAR(math::bernoulli::quadratic_coefficient<0>(squared_argument)[0], sum, 1.0e-7) << squared_argument;
    }
}

TEST(bernoulli, quadratic_coefficient_derivatives) {
    // Central differences of the function and of its first derivative
    for (const double squared_argument :
            {-25.0, -M_PI * M_PI, -4.0, -1.002, -0.998, -0.3, -1.0e-4, 0.0, 0.3, 0.998, 1.002, 4.0, 25.0}) {
        const double step = 1.0e-5;
        const std::array<double, 3> q = math::bernoulli::quadratic_coefficient<2>(squared_argument);
        const std::array<double, 3> q_plus = math::bernoulli::quadratic_coefficient<2>(squared_argument + step);
        const std::array<double, 3> q_minus = math::bernoulli::quadratic_coefficient<2>(squared_argument - step);
        EXPECT_NEAR(q[1], (q_plus[0] - q_minus[0]) / (2.0 * step), 1.0e-9) << squared_argument;
        EXPECT_NEAR(q[2], (q_plus[1] - q_minus[1]) / (2.0 * step), 1.0e-9) << squared_argument;
        // Fewer derivatives are the same values
        EXPECT_EQ(math::bernoulli::quadratic_coefficient<1>(squared_argument)[1], q[1]);
        EXPECT_EQ(math::bernoulli::quadratic_coefficient<0>(squared_argument)[0], q[0]);
    }
}

TEST(bernoulli, quadratic_coefficient_continuous_between_series_and_closed_forms) {
    for (const double squared_argument : {-1.0, 1.0}) {
        const std::array<double, 3> below = math::bernoulli::quadratic_coefficient<2>(squared_argument - 1.0e-12);
        const std::array<double, 3> above = math::bernoulli::quadratic_coefficient<2>(squared_argument + 1.0e-12);
        for (std::size_t i = 0; i < below.size(); ++i) {
            EXPECT_NEAR(below[i], above[i], 1.0e-11);
        }
    }
}
