#ifndef MATHBOX_IMPL_INTEGRATOR_HPP
#define MATHBOX_IMPL_INTEGRATOR_HPP

#include <array>
#include <chrono>
#include <numeric>

#include "mathbox/integrator.hpp"
#include "mathbox/lerp.hpp"

namespace math {

template<typename MathType, typename IndependentVariableType>
inline auto Integrator<MathType, IndependentVariableType>::integrate(const IndependentVariableType start,
        const IndependentVariableType end, const IndependentVariableDifferenceType integration_step,
        const MathType& initial_value,
        const IntegrandFunction<MathType, IndependentVariableType>& integrand) const -> MathType {
    // Error handling and preprocessing
    if (integration_step <= IndependentVariableDifferenceType(0)) {
        throw std::runtime_error("Must integrate with a positive integration_step.");
    }
    const bool forwards = start <= end;
    const IndependentVariableDifferenceType directed_integration_step = forwards ? integration_step : -integration_step;
    // Interval integration bounds
    IndependentVariableType current_start = start;
    IndependentVariableType current_end = forwards ? std::min(current_start + directed_integration_step, end)
                                                   : std::max(current_start + directed_integration_step, end);
    MathType integral = initial_value;
    while ((forwards && current_start < end) || (!forwards && current_start > end)) {
        // Integrate over integration interval
        integral = integrate(current_start, current_end, integral, integrand);
        // Shift the integration interval
        current_start = current_end;
        current_end = forwards ? std::min(current_start + directed_integration_step, end)
                               : std::max(current_start + directed_integration_step, end);
    }
    return integral;
}

namespace newton_cotes {

template<int N, typename MathType, typename IndependentVariableType>
inline auto Integrator<N, MathType, IndependentVariableType>::integrate(const IndependentVariableType start,
        const IndependentVariableType end, const MathType& initial_value,
        const IntegrandFunction<MathType, IndependentVariableType>& integrand) const -> MathType {
    MathType integral = initial_value;
    const ArithmeticTypeScalar step_size_ = step_size(start, end);
    for (std::size_t i = 0; i <= N; ++i) {
        integral += step_size_ * weight(i) * integrand(lerp(start, end, alpha(i)));
    }
    return integral;
}

namespace closed {

template<int N, typename MathType, typename IndependentVariableType>
inline auto Integrator<N, MathType, IndependentVariableType>::alpha(const std::size_t i) const -> ArithmeticTypeScalar {
    return static_cast<ArithmeticTypeScalar>(i) / static_cast<ArithmeticTypeScalar>(N);
}

template<int N, typename MathType, typename IndependentVariableType>
inline auto Integrator<N, MathType, IndependentVariableType>::step_size(const IndependentVariableType start,
        const IndependentVariableType end) const -> ArithmeticTypeScalar {
    return to_sec(end - start) / static_cast<ArithmeticTypeScalar>(N);
}

template<int N, typename MathType, typename IndependentVariableType>
inline auto Integrator<N, MathType, IndependentVariableType>::weight(
        const std::size_t i) const -> ArithmeticTypeScalar {
    return Weights<N, ArithmeticTypeScalar>::weight(i);
}

template<typename Scalar>
auto Weights<1, Scalar>::weight(const std::size_t i) -> Scalar {
    static const std::array<Scalar, 2> weights{{static_cast<Scalar>(1.0 / 2.0), static_cast<Scalar>(1.0 / 2.0)}};
    return weights.at(i);
}

template<typename Scalar>
auto Weights<2, Scalar>::weight(const std::size_t i) -> Scalar {
    static const std::array<Scalar, 3> weights{
            {static_cast<Scalar>(1.0 / 3.0), static_cast<Scalar>(4.0 / 3.0), static_cast<Scalar>(1.0 / 3.0)}};
    return weights.at(i);
}

template<typename Scalar>
auto Weights<3, Scalar>::weight(const std::size_t i) -> Scalar {
    static const std::array<Scalar, 4> weights{{static_cast<Scalar>(3.0 / 8.0), static_cast<Scalar>(9.0 / 8.0),
            static_cast<Scalar>(9.0 / 8.0), static_cast<Scalar>(3.0 / 8.0)}};
    return weights.at(i);
}

template<typename Scalar>
auto Weights<4, Scalar>::weight(const std::size_t i) -> Scalar {
    static const std::array<Scalar, 5> weights{{static_cast<Scalar>(14.0 / 45.0), static_cast<Scalar>(64.0 / 45.0),
            static_cast<Scalar>(24.0 / 45.0), static_cast<Scalar>(64.0 / 45.0), static_cast<Scalar>(14.0 / 45.0)}};
    return weights.at(i);
}

}

namespace open {

template<int N, typename MathType, typename IndependentVariableType>
inline auto Integrator<N, MathType, IndependentVariableType>::alpha(const std::size_t i) const -> ArithmeticTypeScalar {
    return static_cast<ArithmeticTypeScalar>(i + 1) / static_cast<ArithmeticTypeScalar>(N + 2);
}

template<int N, typename MathType, typename IndependentVariableType>
inline auto Integrator<N, MathType, IndependentVariableType>::step_size(const IndependentVariableType start,
        const IndependentVariableType end) const -> ArithmeticTypeScalar {
    return to_sec(end - start) / static_cast<ArithmeticTypeScalar>(N + 2);
}

template<int N, typename MathType, typename IndependentVariableType>
inline auto Integrator<N, MathType, IndependentVariableType>::weight(
        const std::size_t i) const -> ArithmeticTypeScalar {
    return Weights<N, ArithmeticTypeScalar>::weight(i);
}

template<typename Scalar>
auto Weights<0, Scalar>::weight(const std::size_t i) -> Scalar {
    static const std::array<Scalar, 1> weights{{static_cast<Scalar>(2.0)}};
    return weights.at(i);
}

template<typename Scalar>
auto Weights<1, Scalar>::weight(const std::size_t i) -> Scalar {
    static const std::array<Scalar, 2> weights{{static_cast<Scalar>(3.0 / 2.0), static_cast<Scalar>(3.0 / 2.0)}};
    return weights.at(i);
}

template<typename Scalar>
auto Weights<2, Scalar>::weight(const std::size_t i) -> Scalar {
    static const std::array<Scalar, 3> weights{
            {static_cast<Scalar>(8.0 / 3.0), static_cast<Scalar>(-4.0 / 3.0), static_cast<Scalar>(8.0 / 3.0)}};
    return weights.at(i);
}

template<typename Scalar>
auto Weights<3, Scalar>::weight(const std::size_t i) -> Scalar {
    static const std::array<Scalar, 4> weights{{static_cast<Scalar>(55.0 / 24.0), static_cast<Scalar>(5.0 / 24.0),
            static_cast<Scalar>(5.0 / 24.0), static_cast<Scalar>(55.0 / 24.0)}};
    return weights.at(i);
}

}

}

}

#endif
