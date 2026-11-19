#ifndef MATHBOX_BERNOULLI_HPP
#define MATHBOX_BERNOULLI_HPP

#include <array>
#include <type_traits>

namespace math::bernoulli {

/**
 * @brief Bernoulli function \f$ B(x) = \frac{x}{e^x - 1} = \sum_{n=0}^{\infty} B_n \frac{x^n}{n!} \f$, the generating
 * function of the Bernoulli numbers \f$ B_n \f$, where \f$ B(0) = 1 \f$. It is evaluated from its quadratic
 * coefficient (`quadratic_coefficient`), which is accurate near zero.
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
 * in exact rational arithmetic, so the result is the floating-point number nearest to \f$ B_n \f$.
 *
 * @tparam Scalar floating-point type of the result
 * @param n at most 30, the largest for which 64-bit integers hold the rational arithmetic
 * @return Scalar
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
 * It is evaluated with its power series for \f$ |z| < 1 \f$, and otherwise with the closed forms in the half
 * argument, which have no removable singularity at an angle of \f$ \pi \f$.
 *
 * @tparam NumDerivatives number of derivatives (at most 2)
 * @tparam Scalar scalar type, which may be an automatic differentiation type
 * @param squared_argument \f$ z = x^2 \f$, which is negative for an imaginary argument (\f$ -\theta^2 \f$ for an
 * angle)
 * @return std::array<Scalar, NumDerivatives + 1> \f$ q \f$, then \f$ \frac{dq}{dz} \f$ and \f$ \frac{d^2 q}{dz^2}
 * \f$ as requested
 */
template<int NumDerivatives, typename Scalar>
    requires(NumDerivatives >= 0 && NumDerivatives <= 2)
std::array<Scalar, NumDerivatives + 1> quadratic_coefficient(const Scalar& squared_argument);

}

#include "mathbox/impl/bernoulli.hpp"

#endif
