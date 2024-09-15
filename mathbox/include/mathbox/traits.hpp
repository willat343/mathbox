#ifndef MATHBOX_TRAITS_HPP
#define MATHBOX_TRAITS_HPP

#include <Eigen/Core>
#include <Eigen/Dense>

namespace math {

template<typename ArithmeticType>
class ArithmeticTypeTraits;

template<>
class ArithmeticTypeTraits<float> {
public:
    using ArithmeticType = float;
    using Scalar = float;

    /**
     * @brief Return zero for this type.
     *
     * @return constexpr ArithmeticType
     */
    static constexpr ArithmeticType zero();
};

template<>
class ArithmeticTypeTraits<double> {
public:
    using ArithmeticType = double;
    using Scalar = double;

    /**
     * @brief Return zero for this type
     *
     * @return constexpr ArithmeticType
     */
    static constexpr ArithmeticType zero();
};

template<>
class ArithmeticTypeTraits<long double> {
public:
    using ArithmeticType = long double;
    using Scalar = long double;

    /**
     * @brief Return zero for this type
     *
     * @return constexpr ArithmeticType
     */
    static constexpr ArithmeticType zero();
};

template<typename Scalar_, int Rows_, int Cols_, int Options_, int MaxRows_, int MaxCols_>
class ArithmeticTypeTraits<Eigen::Matrix<Scalar_, Rows_, Cols_, Options_, MaxRows_, MaxCols_>> {
public:
    using ArithmeticType = Eigen::Matrix<Scalar_, Rows_, Cols_, Options_, MaxRows_, MaxCols_>;
    using Scalar = Scalar_;

    /**
     * @brief Return zero for this type. Note doxygen bug: https://github.com/doxygen/doxygen/issues/8497
     *
     * @return ArithmeticType
     */
    static constexpr ArithmeticType zero();
};

}

#include "mathbox/impl/traits.hpp"

#endif
