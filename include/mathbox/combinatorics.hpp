#ifndef MATHBOX_COMBINATORICS_HPP
#define MATHBOX_COMBINATORICS_HPP

#include <type_traits>

namespace math {

/**
 * @brief Factorial \f$ n! = 1 \cdot 2 \cdots n \f$, where \f$ 0! = 1 \f$.
 *
 * @tparam Scalar arithmetic type of the result, where 64-bit unsigned integers hold up to \f$ 20! \f$
 * @param n
 * @return Scalar
 */
template<typename Scalar>
    requires(std::is_arithmetic_v<Scalar>)
constexpr Scalar factorial(const unsigned int n);

}

#include "mathbox/impl/combinatorics.hpp"

#endif
