#ifndef MATHBOX_TEST_IMPL_DUAL_HPP
#define MATHBOX_TEST_IMPL_DUAL_HPP

#include <cmath>

#include "test/dual.hpp"

namespace dual {

inline Number::Number(const double value, const double derivative)
    : value(value),
      derivative(derivative) {}

inline Number operator+(const Number& lhs, const Number& rhs) {
    return Number(lhs.value + rhs.value, lhs.derivative + rhs.derivative);
}

inline Number operator-(const Number& lhs, const Number& rhs) {
    return Number(lhs.value - rhs.value, lhs.derivative - rhs.derivative);
}

inline Number operator-(const Number& x) {
    return Number(-x.value, -x.derivative);
}

inline Number operator*(const Number& lhs, const Number& rhs) {
    return Number(lhs.value * rhs.value, lhs.derivative * rhs.value + lhs.value * rhs.derivative);
}

inline Number operator/(const Number& lhs, const Number& rhs) {
    return Number(lhs.value / rhs.value,
            (lhs.derivative * rhs.value - lhs.value * rhs.derivative) / (rhs.value * rhs.value));
}

inline bool operator<(const Number& lhs, const Number& rhs) {
    return lhs.value < rhs.value;
}

inline bool operator>(const Number& lhs, const Number& rhs) {
    return lhs.value > rhs.value;
}

inline Number cos(const Number& x) {
    return Number(std::cos(x.value), -std::sin(x.value) * x.derivative);
}

inline Number cosh(const Number& x) {
    return Number(std::cosh(x.value), std::sinh(x.value) * x.derivative);
}

inline Number exp(const Number& x) {
    return Number(std::exp(x.value), std::exp(x.value) * x.derivative);
}

inline Number sin(const Number& x) {
    return Number(std::sin(x.value), std::cos(x.value) * x.derivative);
}

inline Number sinh(const Number& x) {
    return Number(std::sinh(x.value), std::cosh(x.value) * x.derivative);
}

inline Number sqrt(const Number& x) {
    return Number(std::sqrt(x.value), x.derivative / (2.0 * std::sqrt(x.value)));
}

}

#endif
