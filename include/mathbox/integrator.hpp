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
 * @tparam IntegrandType arithmetic type, supporting addition and multiplication with scalar values
 * @tparam IndependentVariableType type of the independent variable (e.g. double, time point)
 */
template<class IntegrandType, IsIndependentVariableType IndependentVariableType>
using IntegrandFunction = std::function<IntegrandType(const IndependentVariableType)>;

template<class IntegratedType_, IsIndependentVariableType IndependentVariableType_,
        class IntegrandType_ = IntegratedType_>
class Integrator {
public:
    using IntegratedType = IntegratedType_;
    using ArithmeticTypeScalar = MathTypeTraits<IntegratedType>::Scalar;
    using IndependentVariableType = IndependentVariableType_;
    using IndependentVariableDifferenceType = decltype(IndependentVariableType() - IndependentVariableType());
    using IntegrandType = IntegrandType_;

    IntegratedType integrate(const IndependentVariableType start, const IndependentVariableType end,
            const IndependentVariableDifferenceType integration_step, const IntegratedType& initial_value,
            const IntegrandFunction<IntegrandType, IndependentVariableType>& integrand) const;

    virtual IntegratedType integrate(const IndependentVariableType start, const IndependentVariableType end,
            const IntegratedType& initial_value,
            const IntegrandFunction<IntegrandType, IndependentVariableType>& integrand) const = 0;
};

namespace newton_cotes {

template<int N_, class IntegratedType_, IsIndependentVariableType IndependentVariableType_,
        class IntegrandType_ = IntegratedType_>
class Integrator : public math::Integrator<IntegratedType_, IndependentVariableType_, IntegrandType_> {
public:
    using Base = math::Integrator<IntegratedType_, IndependentVariableType_, IntegrandType_>;
    using IntegratedType = Base::IntegratedType;
    using ArithmeticTypeScalar = Base::ArithmeticTypeScalar;
    using IndependentVariableType = Base::IndependentVariableType;
    using IntegrandType = Base::IntegrandType;
    static constexpr int N = N_;

    virtual ArithmeticTypeScalar alpha(const std::size_t i) const = 0;

    using Base::integrate;

    IntegratedType integrate(const IndependentVariableType start, const IndependentVariableType end,
            const IntegratedType& initial_value,
            const IntegrandFunction<IntegrandType_, IndependentVariableType_>& integrand) const override;

    virtual ArithmeticTypeScalar step_size(const IndependentVariableType start,
            const IndependentVariableType end) const = 0;

    virtual ArithmeticTypeScalar weight(const std::size_t i) const = 0;
};

namespace closed {

template<int N_, class IntegratedType_, IsIndependentVariableType IndependentVariableType_,
        class IntegrandType_ = IntegratedType_>
class Integrator
    : public math::newton_cotes::Integrator<N_, IntegratedType_, IndependentVariableType_, IntegrandType_> {
public:
    using Base = math::newton_cotes::Integrator<N_, IntegratedType_, IndependentVariableType_, IntegrandType_>;
    using IntegratedType = Base::IntegratedType;
    using ArithmeticTypeScalar = Base::ArithmeticTypeScalar;
    using IndependentVariableType = Base::IndependentVariableType;
    using IntegrandType = Base::IntegrandType;
    static constexpr int N = Base::N;
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
    static_assert(N > 0, "N must be a positive integer.");
    using Scalar = Scalar_;
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

template<int N_, class IntegratedType_, IsIndependentVariableType IndependentVariableType_,
        class IntegrandType_ = IntegratedType_>
class Integrator
    : public math::newton_cotes::Integrator<N_, IntegratedType_, IndependentVariableType_, IntegrandType_> {
public:
    using Base = math::newton_cotes::Integrator<N_, IntegratedType_, IndependentVariableType_, IntegrandType_>;
    using IntegratedType = Base::IntegratedType;
    using ArithmeticTypeScalar = Base::ArithmeticTypeScalar;
    using IndependentVariableType = Base::IndependentVariableType;
    using IntegrandType = Base::IntegrandType;
    static constexpr int N = Base::N;
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
    static_assert(N >= 0, "N must be an integer >= 0.");
    using Scalar = Scalar_;
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

template<class IntegratedType, typename IndependentVariableType, class IntegrandType = IntegratedType>
using Integrator = newton_cotes::open::Integrator<0, IntegratedType, IndependentVariableType, IntegrandType>;

}

namespace open1 {

template<class IntegratedType, typename IndependentVariableType, class IntegrandType = IntegratedType>
using Integrator = newton_cotes::open::Integrator<1, IntegratedType, IndependentVariableType, IntegrandType>;

}

namespace milnes {

template<class IntegratedType, typename IndependentVariableType, class IntegrandType = IntegratedType>
using Integrator = newton_cotes::open::Integrator<2, IntegratedType, IndependentVariableType, IntegrandType>;

}

namespace open3 {

template<class IntegratedType, typename IndependentVariableType, class IntegrandType = IntegratedType>
using Integrator = newton_cotes::open::Integrator<3, IntegratedType, IndependentVariableType, IntegrandType>;

}

namespace trapezoidal {

template<class IntegratedType, typename IndependentVariableType, class IntegrandType = IntegratedType>
using Integrator = newton_cotes::closed::Integrator<1, IntegratedType, IndependentVariableType, IntegrandType>;

}

namespace simpsons {

template<class IntegratedType, typename IndependentVariableType, class IntegrandType = IntegratedType>
using Integrator = newton_cotes::closed::Integrator<2, IntegratedType, IndependentVariableType, IntegrandType>;

}

namespace simpsons38 {

template<class IntegratedType, typename IndependentVariableType, class IntegrandType = IntegratedType>
using Integrator = newton_cotes::closed::Integrator<3, IntegratedType, IndependentVariableType, IntegrandType>;

}

namespace booles {

template<class IntegratedType, typename IndependentVariableType, class IntegrandType = IntegratedType>
using Integrator = newton_cotes::closed::Integrator<4, IntegratedType, IndependentVariableType, IntegrandType>;

}

}

}

#include "mathbox/impl/integrator.hpp"

#endif
