#ifndef MATHBOX_ROLLING_SHAPE_HPP
#define MATHBOX_ROLLING_SHAPE_HPP

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <concepts>
#include <cppbox/enum.hpp>
#include <cppbox/ranges.hpp>

namespace math {

/**
 * @brief Rolling shape types.
 *
 * SPHERE corresponds to `Sphere`.
 * TORUS corresponds to `Torus`.
 */
CREATE_SMART_ENUM(RollingShapeType, SPHERE, TORUS)

/**
 * @brief Shape of a rigid body that rolls without slipping on a plane, described about its centre, a point fixed in
 * the body.
 *
 * Rolling without slipping holds the material point at the contact at rest, so the centre \f$ \mathbf{c} \f$ moves
 * with the angular velocity \f$ \boldsymbol{\omega} \f$ of the body as
 * \f$ \dot{\mathbf{c}} = \boldsymbol{\omega} \times \boldsymbol{\rho} \f$, where
 * \f$ \boldsymbol{\rho} = \mathbf{c} - \mathbf{p} \f$ is the vector from the contact point \f$ \mathbf{p} \f$ to the
 * centre. Vectors and orientations are expressed in a frame fixed to the plane.
 *
 * Shapes:
 * 1. `Sphere` (SPHERE), e.g. a spherical foot.
 * 2. `Torus` (TORUS), e.g. a wheel with a tyre of circular cross-section, or a thin disc.
 *
 * Use:
 * 1. Each shape computes \f$ \boldsymbol{\rho} \f$ with its `contact_to_centre`, templated on the scalar type of the
 *    evaluation, which may be an automatic differentiation type.
 * 2. `rolling_displacement` integrates the rolling of any shape between two orientations, or through a sequence of
 *    orientations.
 * 3. A shape held as a `RollingShape`, e.g. chosen at run time from a configuration, is recovered with `type()` and
 *    `as_sphere` or `as_torus` before evaluation.
 *
 * A spherical foot that rolls by 0.1 rad about the x axis on level ground moves its centre by 3.3 mm along -y:
 * @code{.cpp}
 * const std::unique_ptr<math::RollingShape<double>> foot = std::make_unique<math::Sphere<double>>(0.033);
 * const Eigen::Vector3d displacement = math::rolling_displacement<double>(math::as_sphere(*foot),
 *         Eigen::Quaterniond::Identity(), Eigen::Quaterniond(Eigen::AngleAxisd(0.1, Eigen::Vector3d::UnitX())),
 *         Eigen::Vector3d::UnitZ());
 * @endcode
 *
 * @tparam Scalar_ floating-point scalar type of the shape's parameters
 */
template<std::floating_point Scalar_ = double>
class RollingShape {
public:
    using Scalar = Scalar_;

    /**
     * @brief Destroy the rolling shape.
     *
     */
    virtual ~RollingShape() = default;

    /**
     * @brief Get the type of the shape.
     *
     * @return RollingShapeType
     */
    virtual RollingShapeType type() const = 0;
};

/**
 * @brief Sphere of radius \f$ r \f$ about its centre, e.g. a spherical foot. On a plane, the centre is above the
 * contact point along the normal, \f$ \boldsymbol{\rho} = r \mathbf{n} \f$.
 *
 * @code{.cpp}
 * const math::Sphere<double> foot(0.033);
 * const Eigen::Vector3d contact_to_centre =
 *         foot.contact_to_centre<double>(Eigen::Quaterniond::Identity(), Eigen::Vector3d::UnitZ());
 * @endcode
 *
 * @tparam Scalar_ floating-point scalar type of the radius
 */
template<std::floating_point Scalar_ = double>
class Sphere final : public RollingShape<Scalar_> {
public:
    /**
     * @brief Construct a sphere.
     *
     * @param radius radius \f$ r \geq 0 \f$
     * @throws std::invalid_argument if the radius is negative
     */
    explicit Sphere(const Scalar_ radius);

    /**
     * @brief Compute the vector from the contact point to the centre, \f$ \boldsymbol{\rho} = r \mathbf{n} \f$.
     *
     * @tparam T scalar type of the evaluation, which may be an automatic differentiation type
     * @param orientation orientation of the body, on which the vector does not depend
     * @param normal unit normal \f$ \mathbf{n} \f$ of the plane, pointing towards the body
     * @return Eigen::Vector<T, 3>
     */
    template<typename T>
    Eigen::Vector<T, 3> contact_to_centre(const Eigen::Quaternion<T>& orientation,
            const Eigen::Vector<T, 3>& normal) const;

    /**
     * @brief Get the radius.
     *
     * @return Scalar_
     */
    Scalar_ radius() const;

    RollingShapeType type() const override;

private:
    // Radius
    Scalar_ radius_;
};

/**
 * @brief Torus about its centre, e.g. a wheel with a tyre of circular cross-section: the circle of radius \f$ R \f$
 * (the major radius) about the unit axis \f$ \mathbf{a} \f$ fixed in the body, swept by a circle of radius \f$ r \f$
 * (the minor radius). A thin disc, e.g. a wheel without a tyre, has \f$ r = 0 \f$.
 *
 * On a plane, the torus touches a distance \f$ r \f$ towards the plane from the point of its major circle closest to
 * the plane, so \f$ \boldsymbol{\rho} = R \mathbf{d} + r \mathbf{n} \f$, where \f$ \mathbf{d} \f$ is the unit
 * direction of the normal within the plane of the major circle.
 *
 * An upright cylinder touches along a line and rolls like a thin disc of the same radius. When tilted, a uniformly
 * stiff flat tread of width \f$ w \f$, deflected by \f$ \delta \f$ under load, moves its centre of pressure along the
 * axis like a torus with \f$ r = w^2 / (12 \delta) \f$ and \f$ R + r \f$ the radius of the wheel, while the whole
 * tread stays in contact.
 *
 * An upright wheel that rolls by 0.5 rad about its axis moves its centre by 60 mm along x:
 * @code{.cpp}
 * const math::Torus<double> wheel(Eigen::Vector3d::UnitY(), 0.1, 0.02);
 * const Eigen::Vector3d displacement = math::rolling_displacement<double>(wheel, Eigen::Quaterniond::Identity(),
 *         Eigen::Quaterniond(Eigen::AngleAxisd(0.5, Eigen::Vector3d::UnitY())), Eigen::Vector3d::UnitZ());
 * @endcode
 *
 * @tparam Scalar_ floating-point scalar type of the axis and radii
 */
template<std::floating_point Scalar_ = double>
class Torus final : public RollingShape<Scalar_> {
public:
    /**
     * @brief Construct a torus.
     *
     * @param axis direction of the axis in the body frame, of any non-zero length
     * @param major_radius major radius \f$ R > 0 \f$
     * @param minor_radius minor radius \f$ r \geq 0 \f$
     * @throws std::invalid_argument if the axis is zero, the major radius is not positive or the minor radius is
     * negative
     */
    Torus(const Eigen::Vector<Scalar_, 3>& axis, const Scalar_ major_radius, const Scalar_ minor_radius);

    /**
     * @brief Get the unit direction of the axis in the body frame.
     *
     * @return const Eigen::Vector<Scalar_, 3>&
     */
    const Eigen::Vector<Scalar_, 3>& axis() const;

    /**
     * @brief Compute the vector from the contact point to the centre,
     * \f$ \boldsymbol{\rho} = R \mathbf{d} + r \mathbf{n} \f$, where \f$ \mathbf{d} \f$ is
     * \f$ \mathbf{n} - (\mathbf{n} \cdot \mathbf{a}) \mathbf{a} \f$ normalised, with the axis \f$ \mathbf{a} \f$ at
     * the orientation of the body.
     *
     * @tparam T scalar type of the evaluation, which may be an automatic differentiation type
     * @param orientation orientation of the body
     * @param normal unit normal \f$ \mathbf{n} \f$ of the plane, pointing towards the body
     * @return Eigen::Vector<T, 3>
     * @throws std::invalid_argument if the axis is parallel to the normal (the torus lies flat on the plane)
     */
    template<typename T>
    Eigen::Vector<T, 3> contact_to_centre(const Eigen::Quaternion<T>& orientation,
            const Eigen::Vector<T, 3>& normal) const;

    /**
     * @brief Get the major radius.
     *
     * @return Scalar_
     */
    Scalar_ major_radius() const;

    /**
     * @brief Get the minor radius.
     *
     * @return Scalar_
     */
    Scalar_ minor_radius() const;

    RollingShapeType type() const override;

private:
    // Unit direction of the axis in the body frame
    Eigen::Vector<Scalar_, 3> axis_;
    // Major radius
    Scalar_ major_radius_;
    // Minor radius
    Scalar_ minor_radius_;
};

/**
 * @brief Shape derived from `RollingShape`.
 *
 * @tparam Shape
 */
template<typename Shape>
concept IsRollingShape = std::derived_from<Shape, RollingShape<typename Shape::Scalar>>;

/**
 * @brief Get a rolling shape as a sphere.
 *
 * @tparam Scalar_
 * @param shape rolling shape of type SPHERE
 * @return const Sphere<Scalar_>&
 * @throws std::bad_cast if the shape is not a sphere
 */
template<std::floating_point Scalar_>
const Sphere<Scalar_>& as_sphere(const RollingShape<Scalar_>& shape);

/**
 * @brief Get a rolling shape as a torus.
 *
 * @tparam Scalar_
 * @param shape rolling shape of type TORUS
 * @return const Torus<Scalar_>&
 * @throws std::bad_cast if the shape is not a torus
 */
template<std::floating_point Scalar_>
const Torus<Scalar_>& as_torus(const RollingShape<Scalar_>& shape);

/**
 * @brief Compute the displacement of the centre of a shape while the body rolls from orientation \f$ \mathbf{q}' \f$
 * to orientation \f$ \mathbf{q}'' \f$, \f$ \int \boldsymbol{\omega} \times \boldsymbol{\rho} \, dt \f$.
 *
 * The integral is approximated by the rotation vector of the increment,
 * \f$ \Delta\boldsymbol{\theta} = \mathrm{Log}(\mathbf{q}'' \mathbf{q}'^{-1}) \f$, crossed with
 * \f$ \boldsymbol{\rho} \f$ at the orientation halfway along the increment,
 * \f$ \mathrm{Exp}(\Delta\boldsymbol{\theta} / 2) \mathbf{q}' \f$. This is exact for a sphere rotating about a fixed
 * axis and for a torus rolling about its axis.
 *
 * @tparam T scalar type of the evaluation, which may be an automatic differentiation type
 * @tparam Shape rolling shape, e.g. `Sphere` or `Torus`
 * @param shape
 * @param start_orientation orientation \f$ \mathbf{q}' \f$ of the body at the start
 * @param end_orientation orientation \f$ \mathbf{q}'' \f$ of the body at the end, less than half a turn from the
 * start
 * @param normal unit normal \f$ \mathbf{n} \f$ of the plane, pointing towards the body
 * @return Eigen::Vector<T, 3>
 */
template<typename T, IsRollingShape Shape>
Eigen::Vector<T, 3> rolling_displacement(const Shape& shape, const Eigen::Quaternion<T>& start_orientation,
        const Eigen::Quaternion<T>& end_orientation, const Eigen::Vector<T, 3>& normal);

/**
 * @brief Compute the displacement of the centre of a shape while the body rolls through the orientations
 * \f$ \mathbf{q}_0, \mathbf{q}_1, \dots, \mathbf{q}_m \f$ in turn, e.g. samples of its orientation over an interval,
 * as the sum of `rolling_displacement` from each orientation to the next.
 *
 * @tparam T scalar type of the evaluation, which may be an automatic differentiation type
 * @tparam Shape rolling shape, e.g. `Sphere` or `Torus`
 * @tparam Orientations forward range of `Eigen::Quaternion<T>`
 * @param shape
 * @param orientations orientations of the body in order, each less than half a turn from the one before
 * @param normal unit normal \f$ \mathbf{n} \f$ of the plane, pointing towards the body
 * @return Eigen::Vector<T, 3> zero for a single orientation
 * @throws std::invalid_argument if there are no orientations
 */
template<typename T, IsRollingShape Shape, cppbox::ForwardRangeOf<Eigen::Quaternion<T>> Orientations>
Eigen::Vector<T, 3> rolling_displacement(const Shape& shape, const Orientations& orientations,
        const Eigen::Vector<T, 3>& normal);

}

#include "mathbox/impl/rolling_shape.hpp"

#endif
