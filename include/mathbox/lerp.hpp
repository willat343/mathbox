#ifndef MATHBOX_LERP_HPP
#define MATHBOX_LERP_HPP

#include "mathbox/traits.hpp"

namespace math {

/**
 * @brief Perform linear interpolation or extrapolation from \f$\mathbf{y}_0\f$ to \f$\mathbf{y}_1\f$, using the
 * formula:
 *
 * \f[
 *      \mathbf{y} = (1 - \alpha) \mathbf{y}_0 + \alpha \mathbf{y}_1
 * \f]
 *
 * If \f$\alpha \in [0, 1]\f$, the function performs interpolation, and otherwise it extrapolates.
 *
 * @tparam ArithmeticType must support multiplication with Scalar, and addition through the + operator
 * @tparam Scalar must be a floating point type
 * @param y_0 start point (\f$\alpha = 0\f$) of interpolation \f$\mathbf{y}_0\f$
 * @param y_1 end point (\f$\alpha = 1\f$) of interpolation \f$\mathbf{y}_1\f$
 * @param alpha interpolation or extrapolation factor \f$\alpha\f$
 * @return Scalar interpolation or extrapolation result \f$\mathbf{y}\f$
 */
template<typename ArithmeticType, typename Scalar = ArithmeticTypeTraits<ArithmeticType>::Scalar>
Scalar lerp(const ArithmeticType& y_0, const ArithmeticType& y_1, const Scalar alpha);

/**
 * @brief Linear function which passes through \f$(x_0, y_0)\f$ and \f$(x_1, y_1)\f$. The functoin returns `lerp(y_0,
 * y_1, alpha)` where `alpha` \f$\alpha\f$ is:
 *
 * \f[
 *      \alpha = \frac{x - x_0}{x_1 - x_0}
 * \f]
 *
 * @tparam ArithmeticType
 * @tparam Scalar
 * @param x
 * @param x_0
 * @param x_1
 * @param y_0
 * @param y_1
 * @return ArithmeticType
 */
template<typename ArithmeticType, typename Scalar = ArithmeticTypeTraits<ArithmeticType>::Scalar>
inline ArithmeticType linear_function(const Scalar x, const Scalar x_0, const Scalar x_1, const ArithmeticType y_0,
        const ArithmeticType y_1);

}

#include "mathbox/impl/lerp.hpp"

#endif
