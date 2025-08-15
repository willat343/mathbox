#ifndef MATHBOX_IMPL_LERP_HPP
#define MATHBOX_IMPL_LERP_HPP

#include <chrono>
#include <cppbox/exceptions.hpp>
#include <stdexcept>

#include "mathbox/lerp.hpp"

namespace math {

template<typename T, typename Scalar>
    requires(is_math_type_v<T> || cppbox::is_time_point_or_duration_v<T>)
inline T lerp(const T& y_0, const T& y_1, const Scalar alpha) {
    if constexpr (is_math_type_v<T>) {
        return (static_cast<Scalar>(1) - alpha) * y_0 + alpha * y_1;
    } else if constexpr (cppbox::is_time_point_v<T>) {
        return y_0 + std::chrono::duration_cast<typename T::duration>(alpha * (y_1 - y_0));
    } else if constexpr (cppbox::is_duration_v<T>) {
        return std::chrono::duration_cast<T>((static_cast<Scalar>(1) - alpha) * y_0 + alpha * y_1);
    }
}

template<IsMathType MathType, IsIndependentVariableType IndependentVariableType>
inline MathType linear_function(const IndependentVariableType x, const IndependentVariableType x_0,
        const IndependentVariableType x_1, const MathType y_0, const MathType y_1) {
    throw_if(x_0 == x_1, "Divide by zero encountered in linear_function because x_0 == x_1.");
    return lerp(y_0, y_1, cppbox::to_sec(x - x_0) / cppbox::to_sec(x_1 - x_0));
}

}

#endif
