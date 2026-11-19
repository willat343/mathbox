#ifndef MATHBOX_IMPL_BERNOULLI_HPP
#define MATHBOX_IMPL_BERNOULLI_HPP

#include <cppbox/exceptions.hpp>
#include <cstdint>
#include <numeric>
#include <stdexcept>

#include "mathbox/bernoulli.hpp"
#include "mathbox/combinatorics.hpp"
#include "mathbox/elementary_functions.hpp"

namespace math::bernoulli {

template<typename Scalar>
Scalar function(const Scalar& x) {
    // B(x) = 1 - x / 2 + x^2 q(x^2)
    const Scalar x_squared = x * x;
    return 1.0 - 0.5 * x + x_squared * quadratic_coefficient<0>(x_squared)[0];
}

template<typename Scalar>
    requires(std::is_floating_point_v<Scalar>)
constexpr Scalar number(const unsigned int n) {
    throw_if(n > 30, "Bernoulli numbers are only available up to B_30.", std::out_of_range);
    // Numerators and denominators of B_0 to B_n, in lowest terms
    std::array<std::int64_t, 31> numerators{1}, denominators{1};
    for (unsigned int m = 1; m <= n; ++m) {
        // sum_{k=0}^{m-1} binomial(m + 1, k) B_k as a fraction in lowest terms
        std::int64_t sum_numerator{0}, sum_denominator{1}, binomial{1};
        for (unsigned int k = 0; k < m; ++k) {
            sum_numerator = sum_numerator * denominators[k] + binomial * numerators[k] * sum_denominator;
            sum_denominator *= denominators[k];
            const std::int64_t divisor = std::gcd(sum_numerator, sum_denominator);
            sum_numerator /= divisor;
            sum_denominator /= divisor;
            // binomial(m + 1, k + 1) = binomial(m + 1, k) (m + 1 - k) / (k + 1)
            binomial = binomial * static_cast<std::int64_t>(m + 1 - k) / static_cast<std::int64_t>(k + 1);
        }
        // B_m = -sum / (m + 1)
        sum_denominator *= static_cast<std::int64_t>(m + 1);
        const std::int64_t divisor = std::gcd(sum_numerator, sum_denominator);
        numerators[m] = -sum_numerator / divisor;
        denominators[m] = sum_denominator / divisor;
    }
    return static_cast<Scalar>(numerators[n]) / static_cast<Scalar>(denominators[n]);
}

template<int NumDerivatives, typename Scalar>
    requires(NumDerivatives >= 0 && NumDerivatives <= 2)
std::array<Scalar, NumDerivatives + 1> quadratic_coefficient(const Scalar& squared_argument) {
    std::array<Scalar, NumDerivatives + 1> q;
    if (squared_argument < 1.0 && squared_argument > -1.0) {
        // Coefficients B_{2k+2} / (2k+2)! of the power series, whose next term is below 1e-16 of the first at 1
        static constexpr std::array<double, 10> series = []() {
            std::array<double, 10> coefficients{};
            for (unsigned int k = 0; k < coefficients.size(); ++k) {
                coefficients[k] = number<double>(2 * k + 2) / factorial<double>(2 * k + 2);
            }
            return coefficients;
        }();
        // Horner's method for the series and its derivatives
        q.fill(Scalar(0.0));
        for (int k = static_cast<int>(series.size()) - 1; k >= 0; --k) {
            q[0] = q[0] * squared_argument + series[k];
            if constexpr (NumDerivatives >= 1) {
                if (k >= 1) {
                    q[1] = q[1] * squared_argument + k * series[k];
                }
            }
            if constexpr (NumDerivatives >= 2) {
                if (k >= 2) {
                    q[2] = q[2] * squared_argument + k * (k - 1) * series[k];
                }
            }
        }
        return q;
    }
    // Closed forms in y, a real argument x (z = x^2, sign = 1) or an angle theta (z = -theta^2, sign = -1), with the
    // hyperbolic sine n and cosine m of y / 2, or the sine and cosine, and the cotangent c = m / n
    const bool is_real = squared_argument > 0.0;
    const double sign = is_real ? 1.0 : -1.0;
    const Scalar y_squared = sign * squared_argument;
    const Scalar y = math::sqrt(y_squared);
    const Scalar n = is_real ? math::sinh(0.5 * y) : math::sin(0.5 * y);
    const Scalar m = is_real ? math::cosh(0.5 * y) : math::cos(0.5 * y);
    const Scalar c = m / n;
    // q = sign (c / (2 y) - 1 / y^2)
    q[0] = sign * (c / (2.0 * y) - 1.0 / y_squared);
    if constexpr (NumDerivatives >= 1) {
        const Scalar n_squared = n * n;
        // dq/dy = sign (2 / y^3 - c / (2 y^2) - 1 / (4 y n^2))
        const Scalar dq_dy = sign * (2.0 / (y * y_squared) - c / (2.0 * y_squared) - 1.0 / (4.0 * y * n_squared));
        // dq/dz = sign dq/dy / (2 y)
        q[1] = sign * dq_dy / (2.0 * y);
        if constexpr (NumDerivatives >= 2) {
            // d2q/dy2 = sign (-6 / y^4 + c / y^3 + 1 / (2 y^2 n^2) + m / (4 y n^3))
            const Scalar d2q_dy2 = sign * (-6.0 / (y_squared * y_squared) + c / (y * y_squared) +
                                                  1.0 / (2.0 * y_squared * n_squared) + m / (4.0 * y * n_squared * n));
            // d2q/dz2 = (d2q/dy2 - dq/dy / y) / (4 y^2)
            q[2] = (d2q_dy2 - dq_dy / y) / (4.0 * y_squared);
        }
    }
    return q;
}

}

#endif
