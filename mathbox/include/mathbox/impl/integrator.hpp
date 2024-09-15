#ifndef MATHBOX_IMPL_INTEGRATOR_HPP
#define MATHBOX_IMPL_INTEGRATOR_HPP

#include <array>
#include <numeric>

#include "mathbox/integrator.hpp"
#include "mathbox/lerp.hpp"

namespace math {

template<typename ArithmeticType, typename Scalar>
inline auto Integrator<ArithmeticType, Scalar>::integrate(const Scalar start, const Scalar end,
        const int num_subintervals,
        const IntegrandFunction<ArithmeticType, Scalar>& integrand) const -> ArithmeticType {
    if (num_subintervals < 1) {
        throw std::runtime_error("Must integrate with a positive number of subintervals, however " +
                                 std::to_string(num_subintervals) + " was given.");
    }
    ArithmeticType integral = ArithmeticTypeTraits<ArithmeticType>::zero();
    for (int i = 0; i < num_subintervals; ++i) {
        integral += integrate(lerp(start, end, static_cast<Scalar>(i) / static_cast<Scalar>(num_subintervals)),
                lerp(start, end, static_cast<Scalar>(i + 1) / static_cast<Scalar>(num_subintervals)), integrand);
    }
    return integral;
}

namespace newton_cotes {

template<int N, typename ArithmeticType, typename Scalar>
inline auto Integrator<N, ArithmeticType, Scalar>::integrate(const Scalar start, const Scalar end,
        const IntegrandFunction<ArithmeticType, Scalar>& integrand) const -> ArithmeticType {
    ArithmeticType integral = ArithmeticTypeTraits<ArithmeticType>::zero();
    for (std::size_t i = 0; i <= N; ++i) {
        integral += weight(i) * integrand(lerp(start, end, alpha(i)));
    }
    return step_size(start, end) * integral;
}

namespace closed {

template<int N, typename ArithmeticType, typename Scalar>
inline auto Integrator<N, ArithmeticType, Scalar>::alpha(const std::size_t i) const -> Scalar {
    return static_cast<Scalar>(i) / static_cast<Scalar>(N);
}

template<int N, typename ArithmeticType, typename Scalar>
inline auto Integrator<N, ArithmeticType, Scalar>::step_size(const Scalar start, const Scalar end) const -> Scalar {
    return (end - start) / static_cast<Scalar>(N);
}

template<int N, typename ArithmeticType, typename Scalar>
inline auto Integrator<N, ArithmeticType, Scalar>::weight(const std::size_t i) const -> Scalar {
    return Weights<N, Scalar>::weight(i);
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

template<int N, typename ArithmeticType, typename Scalar>
inline auto Integrator<N, ArithmeticType, Scalar>::alpha(const std::size_t i) const -> Scalar {
    return static_cast<Scalar>(i + 1) / static_cast<Scalar>(N + 2);
}

template<int N, typename ArithmeticType, typename Scalar>
inline auto Integrator<N, ArithmeticType, Scalar>::step_size(const Scalar start, const Scalar end) const -> Scalar {
    return (end - start) / static_cast<Scalar>(N + 2);
}

template<int N, typename ArithmeticType, typename Scalar>
inline auto Integrator<N, ArithmeticType, Scalar>::weight(const std::size_t i) const -> Scalar {
    return Weights<N, Scalar>::weight(i);
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
