#ifndef MATHBOX_IMPL_ROLLING_SHAPE_HPP
#define MATHBOX_IMPL_ROLLING_SHAPE_HPP

#include <cppbox/exceptions.hpp>
#include <functional>
#include <iterator>
#include <numeric>
#include <ranges>
#include <stdexcept>
#include <vector>

#include "mathbox/geometry.hpp"
#include "mathbox/rolling_shape.hpp"

namespace math {

template<std::floating_point Scalar_>
inline Sphere<Scalar_>::Sphere(const Scalar_ radius)
    : radius_(radius) {
    throw_if(radius < Scalar_(0), "The radius of a sphere must be non-negative.", std::invalid_argument);
}

template<std::floating_point Scalar_>
template<typename T>
Eigen::Vector<T, 3> Sphere<Scalar_>::contact_to_centre(const Eigen::Quaternion<T>&,
        const Eigen::Vector<T, 3>& normal) const {
    // ρ = r n
    return T(radius_) * normal;
}

template<std::floating_point Scalar_>
inline Scalar_ Sphere<Scalar_>::radius() const {
    return radius_;
}

template<std::floating_point Scalar_>
inline RollingShapeType Sphere<Scalar_>::type() const {
    return RollingShapeType::SPHERE;
}

template<std::floating_point Scalar_>
inline Torus<Scalar_>::Torus(const Eigen::Vector<Scalar_, 3>& axis, const Scalar_ major_radius,
        const Scalar_ minor_radius)
    : axis_(axis.normalized()),
      major_radius_(major_radius),
      minor_radius_(minor_radius) {
    throw_if(axis.isZero(), "The axis of a torus must be non-zero.", std::invalid_argument);
    throw_if(major_radius <= Scalar_(0), "The major radius of a torus must be positive.", std::invalid_argument);
    throw_if(minor_radius < Scalar_(0), "The minor radius of a torus must be non-negative.", std::invalid_argument);
}

template<std::floating_point Scalar_>
inline const Eigen::Vector<Scalar_, 3>& Torus<Scalar_>::axis() const {
    return axis_;
}

template<std::floating_point Scalar_>
template<typename T>
Eigen::Vector<T, 3> Torus<Scalar_>::contact_to_centre(const Eigen::Quaternion<T>& orientation,
        const Eigen::Vector<T, 3>& normal) const {
    // Component of the normal within the plane of the major circle, n - (n · a) a, with the axis a at the orientation
    const Eigen::Vector<T, 3> oriented_axis = orientation * axis_.template cast<T>();
    const Eigen::Vector<T, 3> in_plane_normal = normal - normal.dot(oriented_axis) * oriented_axis;
    throw_if(in_plane_normal.isZero(), "The torus lies flat on the plane, so its contact point is not unique.",
            std::invalid_argument);
    // ρ = R d + r n, with d = (n - (n · a) a) / |n - (n · a) a|
    return T(major_radius_) * in_plane_normal.normalized() + T(minor_radius_) * normal;
}

template<std::floating_point Scalar_>
inline Scalar_ Torus<Scalar_>::major_radius() const {
    return major_radius_;
}

template<std::floating_point Scalar_>
inline Scalar_ Torus<Scalar_>::minor_radius() const {
    return minor_radius_;
}

template<std::floating_point Scalar_>
inline RollingShapeType Torus<Scalar_>::type() const {
    return RollingShapeType::TORUS;
}

template<std::floating_point Scalar_>
inline const Sphere<Scalar_>& as_sphere(const RollingShape<Scalar_>& shape) {
    return dynamic_cast<const Sphere<Scalar_>&>(shape);
}

template<std::floating_point Scalar_>
inline const Torus<Scalar_>& as_torus(const RollingShape<Scalar_>& shape) {
    return dynamic_cast<const Torus<Scalar_>&>(shape);
}

template<typename T, IsRollingShape Shape>
Eigen::Vector<T, 3> rolling_displacement(const Shape& shape, const Eigen::Quaternion<T>& start_orientation,
        const Eigen::Quaternion<T>& end_orientation, const Eigen::Vector<T, 3>& normal) {
    // Rotation vector of the increment, Δθ = Log(q'' q'^-1), the shortest rotation
    const Eigen::Vector<T, 3> increment =
            quaternion_to_rotation_vector(Eigen::Quaternion<T>(end_orientation * start_orientation.conjugate()));
    // Δθ × ρ, with ρ at the orientation halfway along the increment, Exp(Δθ / 2) q'
    return increment.cross(shape.template contact_to_centre<T>(
            rotation_vector_to_quaternion<T>(T(0.5) * increment) * start_orientation, normal));
}

template<typename T, IsRollingShape Shape, cppbox::ForwardRangeOf<Eigen::Quaternion<T>> Orientations>
Eigen::Vector<T, 3> rolling_displacement(const Shape& shape, const Orientations& orientations,
        const Eigen::Vector<T, 3>& normal) {
    throw_if(std::ranges::empty(orientations), "The rolling displacement requires at least one orientation.",
            std::invalid_argument);
    // Sum of the displacements from q_{i-1} to q_i, for i = 1, ..., m
    return std::inner_product(std::ranges::begin(orientations),
            std::ranges::next(std::ranges::begin(orientations), std::ranges::distance(orientations) - 1),
            std::ranges::next(std::ranges::begin(orientations)), Eigen::Vector<T, 3>::Zero().eval(), std::plus<>(),
            [&shape, &normal](const Eigen::Quaternion<T>& start_orientation,
                    const Eigen::Quaternion<T>& end_orientation) -> Eigen::Vector<T, 3> {
                return rolling_displacement<T>(shape, start_orientation, end_orientation, normal);
            });
}

}

#if !MATHBOX_HEADER_ONLY
namespace math {

extern template class RollingShape<double>;
extern template class Sphere<double>;
extern template class Torus<double>;

extern template Eigen::Vector3d Sphere<double>::contact_to_centre<double>(const Eigen::Quaterniond&,
        const Eigen::Vector3d&) const;

extern template Eigen::Vector3d Torus<double>::contact_to_centre<double>(const Eigen::Quaterniond&,
        const Eigen::Vector3d&) const;

extern template Eigen::Vector3d rolling_displacement<double, Sphere<double>>(const Sphere<double>&,
        const Eigen::Quaterniond&, const Eigen::Quaterniond&, const Eigen::Vector3d&);
extern template Eigen::Vector3d rolling_displacement<double, Torus<double>>(const Torus<double>&,
        const Eigen::Quaterniond&, const Eigen::Quaterniond&, const Eigen::Vector3d&);

extern template Eigen::Vector3d rolling_displacement<double, Sphere<double>, std::vector<Eigen::Quaterniond>>(
        const Sphere<double>&, const std::vector<Eigen::Quaterniond>&, const Eigen::Vector3d&);
extern template Eigen::Vector3d rolling_displacement<double, Torus<double>, std::vector<Eigen::Quaterniond>>(
        const Torus<double>&, const std::vector<Eigen::Quaterniond>&, const Eigen::Vector3d&);

}
#endif

#endif
