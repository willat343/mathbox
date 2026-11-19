#include "mathbox/bernoulli.hpp"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <stdexcept>

#include "test/dual.hpp"

TEST(bernoulli, function) {
    EXPECT_DOUBLE_EQ(math::bernoulli::function(0.0), 1.0);
    // Arguments near zero, either side of the change from the quadratic coefficient to the quotient, and large
    for (const double x : {-700.0, -5.0, -1.001, -0.999, -0.3, -1.0e-9, 1.0e-9, 0.3, 0.999, 1.001, 4.0, 20.0, 700.0}) {
        const double expected = x / std::expm1(x);
        EXPECT_NEAR(math::bernoulli::function(x), expected, 1.0e-14 * expected) << x;
    }
    // No overflow where the exponential overflows
    EXPECT_EQ(math::bernoulli::function(1000.0), 0.0);
    EXPECT_EQ(math::bernoulli::function(-1000.0), 1000.0);
}

TEST(bernoulli, function_derivative) {
    // dB/dx = B(x) (1 - B(x)) / x - B(x) with an automatic differentiation type, and dB/dx = B_1 at zero
    EXPECT_DOUBLE_EQ(math::bernoulli::function(dual::Number(0.0, 1.0)).derivative, -0.5);
    for (const double x : {-5.0, -0.3, 0.3, 4.0}) {
        const dual::Number B = math::bernoulli::function(dual::Number(x, 1.0));
        EXPECT_NEAR(B.derivative, B.value * (1.0 - B.value) / x - B.value, 1.0e-13) << x;
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
    EXPECT_EQ(math::bernoulli::number<float>(4), -1.0f / 30.0f);
}

TEST(bernoulli, number_from_zeta) {
    // B_2n = (-1)^(n+1) 2 (2n)! zeta(2n) / (2 pi)^(2n), with the zeta function summed from its smallest terms
    for (unsigned int n = 2; n <= 15; ++n) {
        double zeta = 0.0;
        for (int k = 2000; k >= 1; --k) {
            zeta += std::pow(static_cast<double>(k), -2.0 * n);
        }
        const double expected =
                (n % 2 == 0 ? -2.0 : 2.0) * std::tgamma(2.0 * n + 1.0) * zeta / std::pow(2.0 * M_PI, 2.0 * n);
        EXPECT_NEAR(math::bernoulli::number<double>(2 * n), expected, 1.0e-10 * std::abs(expected)) << n;
    }
}

TEST(bernoulli, quadratic_coefficient_values) {
    // Values of the power series in exact rational arithmetic, in the series and in both closed forms
    const std::array<std::array<double, 4>, 4> expected{{
            {-4.0, 0.08947684601641732, -0.0016992357303663318, 9.078316894333242e-05},
            {-0.5, 0.08403614963501425, -0.001422588317696492, 6.868170173101179e-05},
            {0.5, 0.0826470540331997, -0.0013564298697824133, 6.371874128558917e-05},
            {4.0, 0.07825882137483282, -0.0011592795770513052, 4.971819862876677e-05},
    }};
    for (const std::array<double, 4>& values : expected) {
        const std::array<double, 3> q = math::bernoulli::quadratic_coefficient<2>(values[0]);
        EXPECT_NEAR(q[0], values[1], 1.0e-14 * std::abs(values[1])) << values[0];
        EXPECT_NEAR(q[1], values[2], 1.0e-13 * std::abs(values[2])) << values[0];
        EXPECT_NEAR(q[2], values[3], 1.0e-12 * std::abs(values[3])) << values[0];
    }
}

TEST(bernoulli, quadratic_coefficient_automatic_differentiation) {
    // With an automatic differentiation type, the derivative of each output is the next output
    for (const double squared_argument : {-25.0, -1.5, -0.3, 0.0, 0.3, 1.5, 25.0}) {
        const std::array<double, 3> q = math::bernoulli::quadratic_coefficient<2>(squared_argument);
        const std::array<dual::Number, 3> dual_q =
                math::bernoulli::quadratic_coefficient<2>(dual::Number(squared_argument, 1.0));
        EXPECT_EQ(dual_q[0].value, q[0]);
        EXPECT_NEAR(dual_q[0].derivative, q[1], 1.0e-12 * std::abs(q[1])) << squared_argument;
        EXPECT_NEAR(dual_q[1].derivative, q[2], 1.0e-10 * std::abs(q[2])) << squared_argument;
    }
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
