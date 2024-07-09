#ifndef MATHBOX_IMPL_GEOMETRY_HPP
#define MATHBOX_IMPL_GEOMETRY_HPP

#include "mathbox/geometry.hpp"

namespace math {

template<typename Scalar>
inline Eigen::AngleAxis<Scalar> to_angleaxis(const Eigen::Matrix<Scalar, 3, 1>& vector) {
    const Scalar angle = vector.norm();
    const Eigen::Matrix<Scalar, 3, 1> axis =
            (angle == static_cast<Scalar>(0) ? Eigen::Matrix<Scalar, 3, 1>::UnitX() : (vector / angle).eval());
    return Eigen::AngleAxis<Scalar>(angle, axis);
}

template<typename Scalar>
inline Eigen::Matrix<Scalar, 3, 1> to_vector(const Eigen::AngleAxis<Scalar>& angleaxis) {
    return angleaxis.angle() * angleaxis.axis();
}

}

#endif
