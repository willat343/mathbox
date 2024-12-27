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
template<typename MathType, typename IndependentVariableType = typename MathTypeTraits<MathType>::Scalar>
using IntegrandFunction = std::function<MathType(const IndependentVariableType)>;

template<typename MathType_, typename IndependentVariableType_ = typename MathTypeTraits<MathType_>::Scalar>
class Integrator {
public:
    using MathType = MathType_;
    static_assert(is_math_type_v<MathType>, "MathType must be a math type.");
    using ArithmeticTypeScalar = MathTypeTraits<MathType>::Scalar;
    using IndependentVariableType = IndependentVariableType_;
    static_assert(std::is_floating_point_v<IndependentVariableType> || is_time_point_v<IndependentVariableType>,
            "IndependentVariableType must be floating point scalar or std::chrono::time_point<Clock, Duration>.");

    MathType integrate(const IndependentVariableType start, const IndependentVariableType end,
            const int num_subintervals, const IntegrandFunction<MathType, IndependentVariableType>& integrand) const;

    virtual MathType integrate(const IndependentVariableType start, const IndependentVariableType end,
            const IntegrandFunction<MathType, IndependentVariableType>& integrand) const = 0;

protected:
    ArithmeticTypeScalar difference_as_scalar(const IndependentVariableType start,
            const IndependentVariableType end) const;
};

namespace newton_cotes {

template<int N_, typename MathType_, typename IndependentVariableType_ = typename MathTypeTraits<MathType_>::Scalar>
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
            const IntegrandFunction<MathType, IndependentVariableType>& integrand) const override;

    virtual ArithmeticTypeScalar step_size(const IndependentVariableType start,
            const IndependentVariableType end) const = 0;

    virtual ArithmeticTypeScalar weight(const std::size_t i) const = 0;
};

namespace closed {

template<int N_, typename MathType_, typename IndependentVariableType_ = typename MathTypeTraits<MathType_>::Scalar>
class Integrator : public math::newton_cotes::Integrator<N_, MathType_, IndependentVariableType_> {
public:
    using Base = math::newton_cotes::Integrator<N_, MathType_, IndependentVariableType_>;
    using MathType = Base::MathType;
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

template<int N_, typename MathType_, typename IndependentVariableType_ = typename MathTypeTraits<MathType_>::Scalar>
class Integrator : public math::newton_cotes::Integrator<N_, MathType_, IndependentVariableType_> {
public:
    using Base = math::newton_cotes::Integrator<N_, MathType_, IndependentVariableType_>;
    using MathType = Base::MathType;
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

template<typename MathType, typename Scalar = typename MathTypeTraits<MathType>::Scalar>
using Integrator = newton_cotes::open::Integrator<0, MathType, Scalar>;

}

namespace open1 {

template<typename MathType, typename Scalar = typename MathTypeTraits<MathType>::Scalar>
using Integrator = newton_cotes::open::Integrator<1, MathType, Scalar>;

}

namespace milnes {

template<typename MathType, typename Scalar = typename MathTypeTraits<MathType>::Scalar>
using Integrator = newton_cotes::open::Integrator<2, MathType, Scalar>;

}

namespace open3 {

template<typename MathType, typename Scalar = typename MathTypeTraits<MathType>::Scalar>
using Integrator = newton_cotes::open::Integrator<3, MathType, Scalar>;

}

namespace trapezoidal {

template<typename MathType, typename Scalar = typename MathTypeTraits<MathType>::Scalar>
using Integrator = newton_cotes::closed::Integrator<1, MathType, Scalar>;

}

namespace simpsons {

template<typename MathType, typename Scalar = typename MathTypeTraits<MathType>::Scalar>
using Integrator = newton_cotes::closed::Integrator<1, MathType, Scalar>;

}

namespace simpsons38 {

template<typename MathType, typename Scalar = typename MathTypeTraits<MathType>::Scalar>
using Integrator = newton_cotes::closed::Integrator<1, MathType, Scalar>;

}

namespace booles {

template<typename MathType, typename Scalar = typename MathTypeTraits<MathType>::Scalar>
using Integrator = newton_cotes::closed::Integrator<1, MathType, Scalar>;

}

}

}

#include "mathbox/impl/integrator.hpp"

#endif
