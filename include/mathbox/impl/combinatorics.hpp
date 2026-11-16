#ifndef MATHBOX_IMPL_COMBINATORICS_HPP
#define MATHBOX_IMPL_COMBINATORICS_HPP

#include "mathbox/combinatorics.hpp"

namespace math {

template<typename Scalar>
    requires(std::is_arithmetic_v<Scalar>)
constexpr Scalar factorial(const unsigned int n) {
    // Product of the integers from 2 to n, which is empty (one) for n < 2
    Scalar product{1};
    for (unsigned int k = 2; k <= n; ++k) {
        product *= static_cast<Scalar>(k);
    }
    return product;
}

}

#endif
