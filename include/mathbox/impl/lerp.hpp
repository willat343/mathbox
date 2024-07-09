#ifndef MATHBOX_IMPL_LERP_HPP
#define MATHBOX_IMPL_LERP_HPP

#include "mathbox/lerp.hpp"

namespace math {

template<typename ArithmeticType, typename Scalar>
Scalar lerp(const ArithmeticType& x, const ArithmeticType& y, const Scalar alpha) {
    return (static_cast<Scalar>(1) - alpha) * x + alpha * y;
}

template<typename ArithmeticType, typename Scalar>
inline ArithmeticType linear_function(const Scalar x, const Scalar x_0, const Scalar x_1, const ArithmeticType y_0,
        const ArithmeticType y_1) {
    if (x_0 == x_1) {
        throw std::runtime_error("Divide by zero encountered in linear_function because x_0 == x_1.");
    }
    return lerp(y_0, y_1, (x - x_0) / (x_1 - x_0));
}

}

#endif
