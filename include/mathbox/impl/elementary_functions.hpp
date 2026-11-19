#ifndef MATHBOX_IMPL_ELEMENTARY_FUNCTIONS_HPP
#define MATHBOX_IMPL_ELEMENTARY_FUNCTIONS_HPP

#include <cmath>

#include "mathbox/elementary_functions.hpp"

namespace math {

template<typename Scalar>
inline Scalar cos(const Scalar& x) {
    // The unqualified call finds std::cos, or the function of the namespace of the scalar type
    using std::cos;
    return cos(x);
}

template<typename Scalar>
inline Scalar cosh(const Scalar& x) {
    // The unqualified call finds std::cosh, or the function of the namespace of the scalar type
    using std::cosh;
    return cosh(x);
}

template<typename Scalar>
inline Scalar sin(const Scalar& x) {
    // The unqualified call finds std::sin, or the function of the namespace of the scalar type
    using std::sin;
    return sin(x);
}

template<typename Scalar>
inline Scalar sinh(const Scalar& x) {
    // The unqualified call finds std::sinh, or the function of the namespace of the scalar type
    using std::sinh;
    return sinh(x);
}

template<typename Scalar>
inline Scalar sqrt(const Scalar& x) {
    // The unqualified call finds std::sqrt, or the function of the namespace of the scalar type
    using std::sqrt;
    return sqrt(x);
}

}

#endif
