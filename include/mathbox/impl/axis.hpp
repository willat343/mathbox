#ifndef MATHBOX_IMPL_AXIS_HPP
#define MATHBOX_IMPL_AXIS_HPP

#include <algorithm>
#include <cmath>

#include "mathbox/axis.hpp"

namespace math {

template<IsVector Derived>
    requires(Derived::SizeAtCompileTime == 3)
inline SignedAxis SignedAxis::for_vector(const Eigen::MatrixBase<Derived>& vector) {
    Eigen::Index axis;
    vector.cwiseAbs().maxCoeff(&axis);
    const bool positive = vector[axis] >= static_cast<typename Derived::Scalar>(0);
    const SignedAxis positive_signed_axis = positive_for(AxisType(axis));
    return positive ? positive_signed_axis : positive_signed_axis.negative();
}

inline SignedAxis SignedAxis::positive_for(const AxisType axis) {
    return SignedAxis(axis * 2);
}

inline SignedAxis SignedAxis::negative_for(const AxisType axis) {
    return SignedAxis(axis * 2 + 1);
}

inline AxisType SignedAxis::axis() const {
    return AxisType(*this / 2);
}

inline bool SignedAxis::is_positive() const {
    return *this % 2 == 0;
}

inline bool SignedAxis::is_negative() const {
    return !is_positive();
}

template<typename Scalar>
inline Scalar SignedAxis::sign() const {
    return is_positive() ? static_cast<Scalar>(1) : static_cast<Scalar>(-1);
}

inline SignedAxis SignedAxis::positive() const {
    return is_positive() ? *this : SignedAxis(*this - 1);
}

inline SignedAxis SignedAxis::negative() const {
    return is_negative() ? *this : SignedAxis(*this + 1);
}

template<typename Scalar>
inline Eigen::Vector<Scalar, 3> SignedAxis::direction() const {
    return Eigen::Vector<Scalar, 3>::Unit(static_cast<Eigen::Index>(axis())) * sign<Scalar>();
}

template<IsVector Derived>
    requires(Derived::SizeAtCompileTime == 3)
inline typename Derived::Scalar SignedAxis::angle_to(const Eigen::MatrixBase<Derived>& vector) const {
    using Scalar = typename Derived::Scalar;
    return std::acos(std::clamp(direction<Scalar>().dot(vector) / vector.norm(), static_cast<Scalar>(-1),
            static_cast<Scalar>(1)));
}

}

#endif
