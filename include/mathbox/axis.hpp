#ifndef MATHBOX_AXIS_HPP
#define MATHBOX_AXIS_HPP

#include <Eigen/Core>
#include <cppbox/enum.hpp>
#include <mathbox/traits.hpp>

namespace math {

// Axis-aligned unsigned direction (one of the 3 principal axes).
CREATE_SMART_ENUM(AxisType, X, Y, Z)

// Axis-aligned signed direction (+/- along each of the 3 principal axes).
CREATE_SMART_ENUM(SignedAxisType, POSITIVE_X, NEGATIVE_X, POSITIVE_Y, NEGATIVE_Y, POSITIVE_Z, NEGATIVE_Z)

/**
 * @brief One of the 6 axis-aligned signed directions (+/- along each of the 3 principal axes), with axis-aware
 * construction and sign queries layered on top of SignedAxisType.
 */
class SignedAxis : public SignedAxisType {
public:
    using SignedAxisType::SignedAxisType;

    /**
     * @brief Get the signed axis a vector (e.g. acceleration) is most consistent with, i.e. the signed axis with the
     * smallest angle to the vector.
     *
     * @tparam Derived
     * @param vector
     * @return SignedAxis
     */
    template<IsVector Derived>
        requires(Derived::SizeAtCompileTime == 3)
    static SignedAxis for_vector(const Eigen::MatrixBase<Derived>& vector);

    /**
     * @brief Get the positive direction for a given axis (e.g. AxisType::X -> POSITIVE_X).
     *
     * @param axis
     * @return SignedAxis
     */
    static SignedAxis positive_for(const AxisType axis);

    /**
     * @brief Get the negative direction for a given axis (e.g. AxisType::X -> NEGATIVE_X).
     *
     * @param axis
     * @return SignedAxis
     */
    static SignedAxis negative_for(const AxisType axis);

    /**
     * @brief Get the principal axis this direction lies along (e.g. POSITIVE_X or NEGATIVE_X -> AxisType::X).
     *
     * @return AxisType
     */
    AxisType axis() const;

    /**
     * @brief Check whether this is a positive direction (e.g. POSITIVE_X).
     *
     * @return true
     */
    bool is_positive() const;

    /**
     * @brief Check whether this is a negative direction (e.g. NEGATIVE_X).
     *
     * @return true
     */
    bool is_negative() const;

    /**
     * @brief Get the sign of this direction (e.g. POSITIVE_X -> +1, NEGATIVE_X -> -1).
     *
     * @tparam Scalar
     * @return Scalar
     */
    template<typename Scalar>
    Scalar sign() const;

    /**
     * @brief Get the positive counterpart of this direction (e.g. NEGATIVE_X -> POSITIVE_X; a no-op if already
     * positive).
     *
     * @return SignedAxis
     */
    SignedAxis positive() const;

    /**
     * @brief Get the negative counterpart of this direction (e.g. POSITIVE_X -> NEGATIVE_X; a no-op if already
     * negative).
     *
     * @return SignedAxis
     */
    SignedAxis negative() const;

    /**
     * @brief Get the unit vector along this signed direction (e.g. POSITIVE_X -> (1,0,0), NEGATIVE_Z -> (0,0,-1)).
     *
     * @tparam Scalar
     * @return Eigen::Vector<Scalar, 3>
     */
    template<typename Scalar>
    Eigen::Vector<Scalar, 3> direction() const;

    /**
     * @brief Get the angle between a vector (e.g. a mean accelerometer reading) and this signed direction, i.e. how
     * far the vector is from being perfectly aligned with this axis.
     *
     * @tparam Derived
     * @param vector
     * @return Derived::Scalar angle, in radians, in [0, pi]
     */
    template<IsVector Derived>
        requires(Derived::SizeAtCompileTime == 3)
    typename Derived::Scalar angle_to(const Eigen::MatrixBase<Derived>& vector) const;
};

}

#include "mathbox/impl/axis.hpp"

#endif
