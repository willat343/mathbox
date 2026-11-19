#ifndef MATHBOX_TEST_DUAL_HPP
#define MATHBOX_TEST_DUAL_HPP

namespace dual {

/**
 * @brief Dual number \f$ a + b \epsilon \f$ with \f$ \epsilon^2 = 0 \f$, a minimal forward-mode automatic
 * differentiation type: a function evaluated at \f$ x + \epsilon \f$ has its derivative at \f$ x \f$ as its dual
 * part. Its functions are in this namespace, as those of automatic differentiation libraries are in theirs, to test
 * functions that are generic in their scalar type.
 */
struct Number {
    /**
     * @brief Construct a dual number.
     *
     * @param value \f$ a \f$
     * @param derivative \f$ b \f$
     */
    Number(const double value = 0.0, const double derivative = 0.0);

    // Real part a
    double value;
    // Dual part b
    double derivative;
};

Number operator+(const Number& lhs, const Number& rhs);

Number operator-(const Number& lhs, const Number& rhs);

Number operator-(const Number& x);

Number operator*(const Number& lhs, const Number& rhs);

Number operator/(const Number& lhs, const Number& rhs);

bool operator<(const Number& lhs, const Number& rhs);

bool operator>(const Number& lhs, const Number& rhs);

Number cos(const Number& x);

Number cosh(const Number& x);

Number exp(const Number& x);

Number sin(const Number& x);

Number sinh(const Number& x);

Number sqrt(const Number& x);

}

#include "test/impl/dual.hpp"

#endif
