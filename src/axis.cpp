#include "mathbox/axis.hpp"

namespace math {

SignedAxis SignedAxis::positive_for(const AxisType axis) {
    return SignedAxis(axis * 2);
}

SignedAxis SignedAxis::negative_for(const AxisType axis) {
    return SignedAxis(axis * 2 + 1);
}

AxisType SignedAxis::axis() const {
    return AxisType(*this / 2);
}

bool SignedAxis::is_positive() const {
    return *this % 2 == 0;
}

bool SignedAxis::is_negative() const {
    return !is_positive();
}

SignedAxis SignedAxis::positive() const {
    return is_positive() ? *this : SignedAxis(*this - 1);
}

SignedAxis SignedAxis::negative() const {
    return is_negative() ? *this : SignedAxis(*this + 1);
}

}
