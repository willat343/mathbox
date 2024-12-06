#ifndef MATHBOX_IMPL_LERP_HPP
#define MATHBOX_IMPL_LERP_HPP

#include <chrono>
#include <stdexcept>
#include <type_traits>

#include "mathbox/lerp.hpp"
#include "mathbox/time.hpp"
#include "mathbox/traits.hpp"

namespace math {

template<typename T, typename Scalar>
inline T lerp(const T& y_0, const T& y_1, const Scalar alpha) {
    static_assert(is_arithmetic_type_v<T> || is_time_point_v<T> || is_duration_v<T>,
            "Type must be ArithmeticType, std::chrono::time_point<Clock, Duration> or std::chrono::duration<Rep, "
            "Period>.");
    if constexpr (is_arithmetic_type_v<T>) {
        return (static_cast<Scalar>(1) - alpha) * y_0 + alpha * y_1;
    } else if constexpr (is_time_point_v<T>) {
        return y_0 + std::chrono::duration_cast<typename T::duration>(alpha * (y_1 - y_0));
    } else if constexpr (is_duration_v<T>) {
        return std::chrono::duration_cast<T>((static_cast<Scalar>(1) - alpha) * y_0 + alpha * y_1);
    }
}

template<typename ArithmeticType, typename IndependentVariableType>
inline ArithmeticType linear_function(const IndependentVariableType x, const IndependentVariableType x_0,
        const IndependentVariableType x_1, const ArithmeticType y_0, const ArithmeticType y_1) {
    static_assert(is_arithmetic_type_v<ArithmeticType>, "ArithmeticType must be arithmetic type.");
    static_assert(std::is_arithmetic_v<IndependentVariableType> || is_time_point_v<IndependentVariableType>,
            "IndependentVariableType must be Scalar or std::chrono::time_point<Clock, Duration>.");
    if (x_0 == x_1) {
        throw std::runtime_error("Divide by zero encountered in linear_function because x_0 == x_1.");
    }
    if constexpr (std::is_arithmetic_v<IndependentVariableType>) {
        return lerp(y_0, y_1, (x - x_0) / (x_1 - x_0));
    } else if constexpr (is_time_point_v<IndependentVariableType>) {
        // Cast is required to avoid integer division.
        using Scalar = ArithmeticTypeTraits<ArithmeticType>::Scalar;
        return lerp(y_0, y_1, to_sec<Scalar>(x - x_0) / to_sec<Scalar>(x_1 - x_0));
    }
}

}

#endif
