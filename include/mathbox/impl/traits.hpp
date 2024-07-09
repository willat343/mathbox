#ifndef MATHBOX_IMPL_TRAITS_HPP
#define MATHBOX_IMPL_TRAITS_HPP

#include "mathbox/traits.hpp"

namespace math {

inline constexpr auto ArithmeticTypeTraits<float>::zero() -> ArithmeticType {
    return 0.0f;
}

inline constexpr auto ArithmeticTypeTraits<double>::zero() -> ArithmeticType {
    return 0.0;
}

inline constexpr auto ArithmeticTypeTraits<long double>::zero() -> ArithmeticType {
    return 0.0l;
}

template<typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
inline constexpr auto ArithmeticTypeTraits<Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>>::zero()
        -> ArithmeticType {
    return ArithmeticType::Zero();
}

}

#endif
