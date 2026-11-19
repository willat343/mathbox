#include "mathbox/elementary_functions.hpp"

#include <gtest/gtest.h>

#include <cmath>

#include "test/dual.hpp"

TEST(elementary_functions, arithmetic_types) {
    // The functions of the standard library, including for float
    const double x = 0.7;
    EXPECT_EQ(math::cos(x), std::cos(x));
    EXPECT_EQ(math::cosh(x), std::cosh(x));
    EXPECT_EQ(math::exp(x), std::exp(x));
    EXPECT_EQ(math::sin(x), std::sin(x));
    EXPECT_EQ(math::sinh(x), std::sinh(x));
    EXPECT_EQ(math::sqrt(x), std::sqrt(x));
    EXPECT_EQ(math::sqrt(0.7f), std::sqrt(0.7f));
}

TEST(elementary_functions, type_with_functions_in_its_namespace) {
    // The functions of the namespace of the scalar type, here with the derivative of each function as the dual part
    const double x = 0.7;
    const dual::Number dual_x(x, 1.0);
    EXPECT_EQ(math::cos(dual_x).value, std::cos(x));
    EXPECT_EQ(math::cos(dual_x).derivative, -std::sin(x));
    EXPECT_EQ(math::cosh(dual_x).derivative, std::sinh(x));
    EXPECT_EQ(math::exp(dual_x).derivative, std::exp(x));
    EXPECT_EQ(math::sin(dual_x).derivative, std::cos(x));
    EXPECT_EQ(math::sinh(dual_x).derivative, std::cosh(x));
    EXPECT_EQ(math::sqrt(dual_x).derivative, 0.5 / std::sqrt(x));
}
