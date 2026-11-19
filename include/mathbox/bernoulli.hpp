#ifndef MATHBOX_BERNOULLI_HPP
#define MATHBOX_BERNOULLI_HPP

#include <array>
#include <type_traits>

namespace math::bernoulli {

/**
 * @brief Bernoulli function \f$ B(x) = \frac{x}{e^x - 1} = \sum_{n=0}^{\infty} B_n \frac{x^n}{n!} \f$, the generating
 * function of the Bernoulli numbers \f$ B_n \f$, where \f$ B(0) = 1 \f$. For \f$ |x| < 1 \f$ it is evaluated from
 * its quadratic coefficient (`quadratic_coefficient`), which is accurate and differentiable at zero, where the
 * quotient is \f$ 0 / 0 \f$.
 *
 * @tparam Scalar scalar type, which may be an automatic differentiation type
 * @param x
 * @return Scalar
 */
template<typename Scalar>
Scalar function(const Scalar& x);

/**
 * @brief Bernoulli number \f$ B_n \f$ (with \f$ B_1 = -\frac{1}{2} \f$), the coefficients of the Bernoulli function
 * \f$ \frac{x}{e^x - 1} = \sum_{n=0}^{\infty} B_n \frac{x^n}{n!} \f$, from the recurrence
 *
 * \f[
 *      B_m = -\frac{1}{m + 1} \sum_{k=0}^{m-1} \binom{m + 1}{k} B_k, \quad B_0 = 1
 * \f]
 *
 * in exact rational arithmetic, so the only rounding is in the conversion of the numerator and denominator of \f$ B_n
 * \f$ to the floating-point type and their division (for `double`, the result is the nearest number to \f$ B_n \f$).
 * The odd Bernoulli numbers after \f$ B_1 \f$ are zero.
 *
 * @tparam Scalar floating-point type of the result
 * @param n at most 30 (the rational arithmetic overflows 64-bit integers for \f$ B_{32} \f$)
 * @return Scalar
 * @throws std::out_of_range if `n` is greater than 30
 */
template<typename Scalar>
    requires(std::is_floating_point_v<Scalar>)
constexpr Scalar number(const unsigned int n);

/**
 * @brief Quadratic coefficient \f$ q \f$ of the Bernoulli function \f$ B(x) = \frac{x}{e^x - 1} \f$ as a function of
 * its squared argument, and its derivatives:
 *
 * \f[
 *      B(x) = 1 - \frac{x}{2} + x^2 q(x^2), \quad q(z) = \sum_{k=0}^{\infty} \frac{B_{2k+2}}{(2k+2)!} z^k =
 *      \sum_{n=1}^{\infty} \frac{2}{4 \pi^2 n^2 + z}
 * \f]
 *
 * with the Bernoulli numbers \f$ B_n \f$. It is smooth at zero, where \f$ q = \frac{1}{12} \f$, and singular at
 * \f$ z = -(2 \pi n)^2 \f$. The squared argument is positive for a real argument \f$ x \f$ and negative for an
 * imaginary argument \f$ x = i \theta \f$, such as an angle \f$ \theta \f$, with the closed forms:
 *
 * \f[
 *      q(x^2) = \frac{1}{x^2} \left(\frac{x}{2} \coth\frac{x}{2} - 1\right), \quad
 *      q(-\theta^2) = \frac{1}{\theta^2} \left(1 - \frac{\theta}{2} \cot\frac{\theta}{2}\right)
 * \f]
 *
 * It gives the Bernoulli function of an operator \f$ \mathbf{A} \f$ with \f$ \mathbf{A}^3 = z \mathbf{A} \f$ as
 * \f$ B(\mathbf{A}) = \mathbf{I} - \frac{1}{2} \mathbf{A} + q(z) \mathbf{A}^2 \f$. An example is the inverse of the
 * Jacobian of the exponential map of a Lie group, \f$ B(\mathrm{ad}) \f$ (up to the sign of the linear term for the
 * left or right Jacobian), of rotations in three dimensions by an angle \f$ \theta \f$, with \f$ z = -\theta^2 \f$.
 * For a real argument, \f$ q(x^2) = L(x / 2) / (2 x) \f$ with the Langevin function \f$ L(x) = \coth x - \frac{1}{x}
 * \f$.
 *
 * It is evaluated with ten terms of its power series for \f$ |z| < 1 \f$, and otherwise with the closed forms above.
 * They are in the half argument because the form of the angle with \f$ \sin\theta \f$ in its denominator is \f$ 0 /
 * 0 \f$ at an angle of \f$ \pi \f$. The relative accuracy of \f$ q \f$ and of its first, second and third derivatives
 * is lowest near \f$ |z| = 1 \f$, at about \f$ 10^{-15} \f$, \f$ 10^{-13} \f$, \f$ 10^{-11} \f$ and \f$ 10^{-9} \f$
 * (from the truncation of the series below, and from cancellation in the closed forms above), and near machine
 * precision away from it.
 *
 * @tparam NumDerivatives number of derivatives (at most 3)
 * @tparam Scalar scalar type, which may be an automatic differentiation type
 * @param squared_argument \f$ z = x^2 \f$, which is negative for an imaginary argument (\f$ -\theta^2 \f$ for an
 * angle)
 * @return std::array<Scalar, NumDerivatives + 1> \f$ q \f$, then \f$ \frac{dq}{dz} \f$, \f$ \frac{d^2 q}{dz^2} \f$
 * and \f$ \frac{d^3 q}{dz^3} \f$ as requested
 */
template<int NumDerivatives, typename Scalar>
    requires(NumDerivatives >= 0 && NumDerivatives <= 3)
std::array<Scalar, NumDerivatives + 1> quadratic_coefficient(const Scalar& squared_argument);

}

#include "mathbox/impl/bernoulli.hpp"

#endif
