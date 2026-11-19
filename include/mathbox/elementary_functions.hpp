#ifndef MATHBOX_ELEMENTARY_FUNCTIONS_HPP
#define MATHBOX_ELEMENTARY_FUNCTIONS_HPP

namespace math {

/**
 * @brief Cosine of a scalar of any type, for functions that are generic in their scalar type. It is `std::cos` for
 * arithmetic types, and otherwise the `cos` function found in the namespace of the scalar type (argument-dependent
 * lookup), e.g. that of an automatic differentiation type.
 *
 * @tparam Scalar
 * @param x
 * @return Scalar
 */
template<typename Scalar>
Scalar cos(const Scalar& x);

/**
 * @brief Hyperbolic cosine of a scalar of any type (see `math::cos`).
 *
 * @tparam Scalar
 * @param x
 * @return Scalar
 */
template<typename Scalar>
Scalar cosh(const Scalar& x);

/**
 * @brief Exponential of a scalar of any type (see `math::cos`).
 *
 * @tparam Scalar
 * @param x
 * @return Scalar
 */
template<typename Scalar>
Scalar exp(const Scalar& x);

/**
 * @brief Sine of a scalar of any type (see `math::cos`).
 *
 * @tparam Scalar
 * @param x
 * @return Scalar
 */
template<typename Scalar>
Scalar sin(const Scalar& x);

/**
 * @brief Hyperbolic sine of a scalar of any type (see `math::cos`).
 *
 * @tparam Scalar
 * @param x
 * @return Scalar
 */
template<typename Scalar>
Scalar sinh(const Scalar& x);

/**
 * @brief Square root of a scalar of any type (see `math::cos`).
 *
 * @tparam Scalar
 * @param x
 * @return Scalar
 */
template<typename Scalar>
Scalar sqrt(const Scalar& x);

}

#include "mathbox/impl/elementary_functions.hpp"

#endif
