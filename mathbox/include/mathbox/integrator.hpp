#ifndef MATHBOX_INTEGRATOR_HPP
#define MATHBOX_INTEGRATOR_HPP

#include <functional>
#include <type_traits>

#include "mathbox/time.hpp"
#include "mathbox/traits.hpp"

namespace math {

/**
 * @brief Integrand function type. std::function is used to allow for a wide variety of callables (e.g. function
 * pointers, lambdas, bound functions, etc.), while still restricting the type.
 *
 * @tparam ArithmeticType arithmetic type, supporting addition and multiplication with scalar values
 * @tparam IndependentVariableType type of the independent variable (e.g. time)
 */
template<typename ArithmeticType,
        typename IndependentVariableType = typename ArithmeticTypeTraits<ArithmeticType>::Scalar>
using IntegrandFunction = std::function<ArithmeticType(const IndependentVariableType)>;

template<typename ArithmeticType_,
        typename IndependentVariableType_ = typename ArithmeticTypeTraits<ArithmeticType_>::Scalar>
class Integrator {
public:
    using ArithmeticType = ArithmeticType_;
    static_assert(is_arithmetic_type_v<ArithmeticType>, "ArithmeticType must be an arithmetic type.");
    using ArithmeticTypeScalar = ArithmeticTypeTraits<ArithmeticType>::Scalar;
    using IndependentVariableType = IndependentVariableType_;
    static_assert(std::is_arithmetic_v<IndependentVariableType> || is_time_point_v<IndependentVariableType>,
            "IndependentVariableType must be an Scalar or std::chrono::time_point<Clock, Duration>.");

    ArithmeticType integrate(const IndependentVariableType start, const IndependentVariableType end,
            const int num_subintervals,
            const IntegrandFunction<ArithmeticType, IndependentVariableType>& integrand) const;

    virtual ArithmeticType integrate(const IndependentVariableType start, const IndependentVariableType end,
            const IntegrandFunction<ArithmeticType, IndependentVariableType>& integrand) const = 0;

protected:
    ArithmeticTypeScalar difference_as_scalar(const IndependentVariableType start,
            const IndependentVariableType end) const;
};

namespace newton_cotes {

template<int N_, typename ArithmeticType_,
        typename IndependentVariableType_ = typename ArithmeticTypeTraits<ArithmeticType_>::Scalar>
class Integrator : public math::Integrator<ArithmeticType_, IndependentVariableType_> {
public:
    using Base = math::Integrator<ArithmeticType_, IndependentVariableType_>;
    using ArithmeticType = Base::ArithmeticType;
    using ArithmeticTypeScalar = Base::ArithmeticTypeScalar;
    using IndependentVariableType = Base::IndependentVariableType;
    static constexpr int N = N_;

    virtual ArithmeticTypeScalar alpha(const std::size_t i) const = 0;

    using Base::integrate;

    ArithmeticType integrate(const IndependentVariableType start, const IndependentVariableType end,
            const IntegrandFunction<ArithmeticType, IndependentVariableType>& integrand) const override;

    virtual ArithmeticTypeScalar step_size(const IndependentVariableType start,
            const IndependentVariableType end) const = 0;

    virtual ArithmeticTypeScalar weight(const std::size_t i) const = 0;
};

namespace closed {

template<int N_, typename ArithmeticType_,
        typename IndependentVariableType_ = typename ArithmeticTypeTraits<ArithmeticType_>::Scalar>
class Integrator : public math::newton_cotes::Integrator<N_, ArithmeticType_, IndependentVariableType_> {
public:
    using Base = math::newton_cotes::Integrator<N_, ArithmeticType_, IndependentVariableType_>;
    using ArithmeticType = Base::ArithmeticType;
    using ArithmeticTypeScalar = Base::ArithmeticTypeScalar;
    using IndependentVariableType = Base::IndependentVariableType;
    static constexpr int N = Base::N;

    // Type Requirements
    static_assert(N > 0, "N must be > 0.");

    ArithmeticTypeScalar alpha(const std::size_t i) const override;

    ArithmeticTypeScalar step_size(const IndependentVariableType start,
            const IndependentVariableType end) const override;

    ArithmeticTypeScalar weight(const std::size_t i) const override;
};

template<int N_, typename Scalar_>
class Weights {
public:
    static constexpr int N = N_;
    using Scalar = Scalar_;

    // Type Requirements
    static_assert(N > 0, "N must be a positive integer.");
    static_assert(std::is_floating_point_v<Scalar>, "Scalar must be floating point type.");

    static Scalar weight(const std::size_t i) = delete;
};

template<typename Scalar_>
class Weights<1, Scalar_> {
public:
    using Scalar = Scalar_;

    static Scalar weight(const std::size_t i);
};

template<typename Scalar_>
class Weights<2, Scalar_> {
public:
    using Scalar = Scalar_;

    static Scalar weight(const std::size_t i);
};

template<typename Scalar_>
class Weights<3, Scalar_> {
public:
    using Scalar = Scalar_;

    static Scalar weight(const std::size_t i);
};

template<typename Scalar_>
class Weights<4, Scalar_> {
public:
    using Scalar = Scalar_;

    static Scalar weight(const std::size_t i);
};

}

namespace open {

template<int N_, typename ArithmeticType_,
        typename IndependentVariableType_ = typename ArithmeticTypeTraits<ArithmeticType_>::Scalar>
class Integrator : public math::newton_cotes::Integrator<N_, ArithmeticType_, IndependentVariableType_> {
public:
    using Base = math::newton_cotes::Integrator<N_, ArithmeticType_, IndependentVariableType_>;
    using ArithmeticType = Base::ArithmeticType;
    using ArithmeticTypeScalar = Base::ArithmeticTypeScalar;
    using IndependentVariableType = Base::IndependentVariableType;
    static constexpr int N = Base::N;

    // Type Requirements
    static_assert(N >= 0, "N must be >= 0.");

    ArithmeticTypeScalar alpha(const std::size_t i) const override;

    ArithmeticTypeScalar step_size(const IndependentVariableType start,
            const IndependentVariableType end) const override;

    ArithmeticTypeScalar weight(const std::size_t i) const override;
};

template<int N_, typename Scalar_>
class Weights {
public:
    static constexpr int N = N_;
    using Scalar = Scalar_;

    // Type Requirements
    static_assert(N >= 0, "N must be an integer >= 0.");
    static_assert(std::is_floating_point_v<Scalar>, "Scalar must be floating point type.");

    static Scalar weight(const std::size_t i) = delete;
};

template<typename Scalar_>
class Weights<0, Scalar_> {
public:
    using Scalar = Scalar_;

    static Scalar weight(const std::size_t i);
};

template<typename Scalar_>
class Weights<1, Scalar_> {
public:
    using Scalar = Scalar_;

    static Scalar weight(const std::size_t i);
};

template<typename Scalar_>
class Weights<2, Scalar_> {
public:
    using Scalar = Scalar_;

    static Scalar weight(const std::size_t i);
};

template<typename Scalar_>
class Weights<3, Scalar_> {
public:
    using Scalar = Scalar_;

    static Scalar weight(const std::size_t i);
};

}

namespace rectangle {

template<typename ArithmeticType, typename Scalar = typename ArithmeticTypeTraits<ArithmeticType>::Scalar>
using Integrator = newton_cotes::open::Integrator<0, ArithmeticType, Scalar>;

}

namespace open1 {

template<typename ArithmeticType, typename Scalar = typename ArithmeticTypeTraits<ArithmeticType>::Scalar>
using Integrator = newton_cotes::open::Integrator<1, ArithmeticType, Scalar>;

}

namespace milnes {

template<typename ArithmeticType, typename Scalar = typename ArithmeticTypeTraits<ArithmeticType>::Scalar>
using Integrator = newton_cotes::open::Integrator<2, ArithmeticType, Scalar>;

}

namespace open3 {

template<typename ArithmeticType, typename Scalar = typename ArithmeticTypeTraits<ArithmeticType>::Scalar>
using Integrator = newton_cotes::open::Integrator<3, ArithmeticType, Scalar>;

}

namespace trapezoidal {

template<typename ArithmeticType, typename Scalar = typename ArithmeticTypeTraits<ArithmeticType>::Scalar>
using Integrator = newton_cotes::closed::Integrator<1, ArithmeticType, Scalar>;

}

namespace simpsons {

template<typename ArithmeticType, typename Scalar = typename ArithmeticTypeTraits<ArithmeticType>::Scalar>
using Integrator = newton_cotes::closed::Integrator<1, ArithmeticType, Scalar>;

}

namespace simpsons38 {

template<typename ArithmeticType, typename Scalar = typename ArithmeticTypeTraits<ArithmeticType>::Scalar>
using Integrator = newton_cotes::closed::Integrator<1, ArithmeticType, Scalar>;

}

namespace booles {

template<typename ArithmeticType, typename Scalar = typename ArithmeticTypeTraits<ArithmeticType>::Scalar>
using Integrator = newton_cotes::closed::Integrator<1, ArithmeticType, Scalar>;

}

}

}

#include "mathbox/impl/integrator.hpp"

#endif
