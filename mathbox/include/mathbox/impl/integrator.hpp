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
        const IndependentVariableType end, const int num_subintervals,
        const IntegrandFunction<MathType, IndependentVariableType>& integrand) const -> MathType {
    if (num_subintervals < 1) {
        throw std::runtime_error("Must integrate with a positive number of subintervals, however " +
                                 std::to_string(num_subintervals) + " was given.");
    }
    MathType integral = MathTypeTraits<MathType>::zero();
    for (int i = 0; i < num_subintervals; ++i) {
        integral += integrate(
                lerp(start, end,
                        static_cast<ArithmeticTypeScalar>(i) / static_cast<ArithmeticTypeScalar>(num_subintervals)),
                lerp(start, end,
                        static_cast<ArithmeticTypeScalar>(i + 1) / static_cast<ArithmeticTypeScalar>(num_subintervals)),
                integrand);
    }
    return integral;
}

template<typename MathType, typename IndependentVariableType>
inline auto Integrator<MathType, IndependentVariableType>::difference_as_scalar(const IndependentVariableType start,
        const IndependentVariableType end) const -> ArithmeticTypeScalar {
    if constexpr (is_time_point_v<IndependentVariableType>) {
        // Convert step size to seconds as a ArithmeticTypeScalar, because integrands are assumed to be per second.
        return to_sec(end - start);
    } else {
        return end - start;
    }
}

namespace newton_cotes {

template<int N, typename MathType, typename IndependentVariableType>
inline auto Integrator<N, MathType, IndependentVariableType>::integrate(const IndependentVariableType start,
        const IndependentVariableType end,
        const IntegrandFunction<MathType, IndependentVariableType>& integrand) const -> MathType {
    MathType integral = MathTypeTraits<MathType>::zero();
    for (std::size_t i = 0; i <= N; ++i) {
        integral += weight(i) * integrand(lerp(start, end, alpha(i)));
    }
    return step_size(start, end) * integral;
}

namespace closed {

template<int N, typename MathType, typename IndependentVariableType>
inline auto Integrator<N, MathType, IndependentVariableType>::alpha(const std::size_t i) const -> ArithmeticTypeScalar {
    return static_cast<ArithmeticTypeScalar>(i) / static_cast<ArithmeticTypeScalar>(N);
}

template<int N, typename MathType, typename IndependentVariableType>
inline auto Integrator<N, MathType, IndependentVariableType>::step_size(const IndependentVariableType start,
        const IndependentVariableType end) const -> ArithmeticTypeScalar {
    return this->difference_as_scalar(start, end) / static_cast<ArithmeticTypeScalar>(N);
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
    return this->difference_as_scalar(start, end) / static_cast<ArithmeticTypeScalar>(N + 2);
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
