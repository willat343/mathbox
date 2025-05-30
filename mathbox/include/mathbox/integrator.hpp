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
 * @tparam MathType arithmetic type, supporting addition and multiplication with scalar values
 * @tparam IndependentVariableType type of the independent variable (e.g. time)
 */
template<IsMathType MathType,
        IsIndependentVariableType IndependentVariableType = typename MathTypeTraits<MathType>::Scalar>
using IntegrandFunction = std::function<MathType(const IndependentVariableType)>;

template<IsMathType MathType_,
        IsIndependentVariableType IndependentVariableType_ = typename MathTypeTraits<MathType_>::Scalar>
class Integrator {
public:
    using MathType = MathType_;
    using ArithmeticTypeScalar = MathTypeTraits<MathType>::Scalar;
    using IndependentVariableType = IndependentVariableType_;
    using IndependentVariableDifferenceType = decltype(IndependentVariableType() - IndependentVariableType());

    MathType integrate(const IndependentVariableType start, const IndependentVariableType end,
            const IndependentVariableDifferenceType integration_step, const MathType& initial_value,
            const IntegrandFunction<MathType, IndependentVariableType>& integrand) const;

    virtual MathType integrate(const IndependentVariableType start, const IndependentVariableType end,
            const MathType& initial_value,
            const IntegrandFunction<MathType, IndependentVariableType>& integrand) const = 0;
};

namespace newton_cotes {

template<int N_, IsMathType MathType_,
        IsIndependentVariableType IndependentVariableType_ = typename MathTypeTraits<MathType_>::Scalar>
class Integrator : public math::Integrator<MathType_, IndependentVariableType_> {
public:
    using Base = math::Integrator<MathType_, IndependentVariableType_>;
    using MathType = Base::MathType;
    using ArithmeticTypeScalar = Base::ArithmeticTypeScalar;
    using IndependentVariableType = Base::IndependentVariableType;
    static constexpr int N = N_;

    virtual ArithmeticTypeScalar alpha(const std::size_t i) const = 0;

    using Base::integrate;

    MathType integrate(const IndependentVariableType start, const IndependentVariableType end,
            const MathType& initial_value,
            const IntegrandFunction<MathType, IndependentVariableType>& integrand) const override;

    virtual ArithmeticTypeScalar step_size(const IndependentVariableType start,
            const IndependentVariableType end) const = 0;

    virtual ArithmeticTypeScalar weight(const std::size_t i) const = 0;
};

namespace closed {

template<int N_, IsMathType MathType_,
        IsIndependentVariableType IndependentVariableType_ = typename MathTypeTraits<MathType_>::Scalar>
class Integrator : public math::newton_cotes::Integrator<N_, MathType_, IndependentVariableType_> {
public:
    using Base = math::newton_cotes::Integrator<N_, MathType_, IndependentVariableType_>;
    using MathType = Base::MathType;
    using ArithmeticTypeScalar = Base::ArithmeticTypeScalar;
    using IndependentVariableType = Base::IndependentVariableType;
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

template<int N_, IsMathType MathType_,
        IsIndependentVariableType IndependentVariableType_ = typename MathTypeTraits<MathType_>::Scalar>
class Integrator : public math::newton_cotes::Integrator<N_, MathType_, IndependentVariableType_> {
public:
    using Base = math::newton_cotes::Integrator<N_, MathType_, IndependentVariableType_>;
    using MathType = Base::MathType;
    using ArithmeticTypeScalar = Base::ArithmeticTypeScalar;
    using IndependentVariableType = Base::IndependentVariableType;
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

template<typename MathType, typename IndependentVariableType = typename MathTypeTraits<MathType>::Scalar>
using Integrator = newton_cotes::open::Integrator<0, MathType, IndependentVariableType>;

}

namespace open1 {

template<typename MathType, typename IndependentVariableType = typename MathTypeTraits<MathType>::Scalar>
using Integrator = newton_cotes::open::Integrator<1, MathType, IndependentVariableType>;

}

namespace milnes {

template<typename MathType, typename IndependentVariableType = typename MathTypeTraits<MathType>::Scalar>
using Integrator = newton_cotes::open::Integrator<2, MathType, IndependentVariableType>;

}

namespace open3 {

template<typename MathType, typename IndependentVariableType = typename MathTypeTraits<MathType>::Scalar>
using Integrator = newton_cotes::open::Integrator<3, MathType, IndependentVariableType>;

}

namespace trapezoidal {

template<typename MathType, typename IndependentVariableType = typename MathTypeTraits<MathType>::Scalar>
using Integrator = newton_cotes::closed::Integrator<1, MathType, IndependentVariableType>;

}

namespace simpsons {

template<typename MathType, typename IndependentVariableType = typename MathTypeTraits<MathType>::Scalar>
using Integrator = newton_cotes::closed::Integrator<2, MathType, IndependentVariableType>;

}

namespace simpsons38 {

template<typename MathType, typename IndependentVariableType = typename MathTypeTraits<MathType>::Scalar>
using Integrator = newton_cotes::closed::Integrator<3, MathType, IndependentVariableType>;

}

namespace booles {

template<typename MathType, typename IndependentVariableType = typename MathTypeTraits<MathType>::Scalar>
using Integrator = newton_cotes::closed::Integrator<4, MathType, IndependentVariableType>;

}

}

}

#include "mathbox/impl/integrator.hpp"

#endif
