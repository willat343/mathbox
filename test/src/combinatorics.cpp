#include "mathbox/combinatorics.hpp"

#include <gtest/gtest.h>

#include <cstdint>

TEST(combinatorics, factorial) {
    // The empty product 0! = 1 and small factorials at compile time, and the largest factorial held by 64-bit unsigned
    // integers, which double also holds exactly
    static_assert(math::factorial<int>(0) == 1 && math::factorial<int>(1) == 1 && math::factorial<int>(5) == 120);
    EXPECT_EQ(math::factorial<std::uint64_t>(20), 2432902008176640000u);
    EXPECT_EQ(math::factorial<double>(20), 2432902008176640000.0);
}
