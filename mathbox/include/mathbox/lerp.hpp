#ifndef MATHBOX_LERP_HPP
#define MATHBOX_LERP_HPP

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
 * To reduce floating point inaccuracies ("catastrophic cancellation" when the input
 * values are large), this formula is preferred over:
 *
 * \f[
 *      \mathbf{y} = \mathbf{y} = \mathbf{y}_0 + \alpha (\mathbf{y}_1 - \mathbf{y}_0)
 * \f]
 *
 * This second formula is used for the linear interpolation of `std::chrono::time_point<Clock, Duration>` objects.
 *
 * @tparam T This must be a mathbox arithmetic type (e.g., floating-point scalars and Eigen matrices), a
 * `std::chrono::time_point<Clock, Duration>` or std::chrono::duration<Rep, Period>` object.
 * @tparam Scalar
 * @param y_0 start point (\f$\alpha = 0\f$) of interpolation \f$\mathbf{y}_0\f$
 * @param y_1 end point (\f$\alpha = 1\f$) of interpolation \f$\mathbf{y}_1\f$
 * @param alpha interpolation or extrapolation factor \f$\alpha\f$
 * @return T interpolation or extrapolation result \f$\mathbf{y}\f$
 */
template<typename T, typename Scalar>
T lerp(const T& y_0, const T& y_1, const Scalar alpha);

/**
 * @brief Linear function which passes through \f$(x_0, y_0)\f$ and \f$(x_1, y_1)\f$. The function returns `lerp(y_0,
 * y_1, alpha)` where `alpha` \f$\alpha\f$ is:
 *
 * \f[
 *      \alpha = \frac{x - x_0}{x_1 - x_0}
 * \f]
 *
 * @tparam ArithmeticType
 * @tparam IndependentVariableType
 * @param x
 * @param x_0
 * @param x_1
 * @param y_0
 * @param y_1
 * @return ArithmeticType
 */
template<typename ArithmeticType, typename IndependentVariableType>
inline ArithmeticType linear_function(const IndependentVariableType x, const IndependentVariableType x_0,
        const IndependentVariableType x_1, const ArithmeticType y_0, const ArithmeticType y_1);

}

#include "mathbox/impl/lerp.hpp"

#endif
