#include "mathbox/rolling_shape.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <stdexcept>
#include <typeinfo>
#include <vector>

/**
 * @brief Check that the rolling displacements of a shape over consecutive 20 ms increments sum to the path integral
 * of \f$ \boldsymbol{\omega} \times \boldsymbol{\rho} \f$, integrated finely while the body rolls with an angular
 * velocity whose axis turns, and that the rolling displacement through the orientations at the ends of the increments
 * is that sum.
 *
 * @tparam Shape
 * @param shape
 */
template<typename Shape>
void check_increments_follow_the_rolling_path(const Shape& shape) {
    const Eigen::Vector3d normal = Eigen::Vector3d(0.1, 0.2, 1.0).normalized();
    const auto angular_velocity = [](const double t) {
        return Eigen::Vector3d(0.5 * std::cos(3.0 * t), 0.8, 0.3 * std::sin(2.0 * t));
    };
    const auto rotation = [](const Eigen::Vector3d& rotation_vector) {
        return Eigen::Quaterniond(Eigen::AngleAxisd(rotation_vector.norm(), rotation_vector.normalized()));
    };
    const double duration = 0.4, fine_step = 1.0e-5, increment = 0.02;
    const int steps_per_increment = static_cast<int>(std::round(increment / fine_step));
    Eigen::Quaterniond orientation = rotation(Eigen::Vector3d(0.1, -0.2, 0.3)), increment_start = orientation;
    std::vector<Eigen::Quaterniond> increment_ends{orientation};
    Eigen::Vector3d path_displacement = Eigen::Vector3d::Zero(), summed_displacement = Eigen::Vector3d::Zero();
    for (int i = 0; i < static_cast<int>(std::round(duration / fine_step)); ++i) {
        // Midpoint rule over each fine step
        const Eigen::Vector3d omega = angular_velocity((i + 0.5) * fine_step);
        path_displacement += fine_step * omega.cross(shape.template contact_to_centre<double>(
                                                 rotation(0.5 * fine_step * omega) * orientation, normal));
        orientation = (rotation(fine_step * omega) * orientation).normalized();
        if ((i + 1) % steps_per_increment == 0) {
            summed_displacement += math::rolling_displacement<double>(shape, increment_start, orientation, normal);
            increment_start = orientation;
            increment_ends.push_back(orientation);
        }
    }
    EXPECT_LT((summed_displacement - path_displacement).norm(), 1.0e-6);
    EXPECT_GT(path_displacement.norm(), 1.0e-2);
    EXPECT_LT((math::rolling_displacement<double>(shape, increment_ends, normal) - summed_displacement).norm(),
            1.0e-15);
}

TEST(rolling_shape, documentation_example) {
    const std::unique_ptr<math::RollingShape<double>> foot = std::make_unique<math::Sphere<double>>(0.033);
    const Eigen::Vector3d displacement =
            math::rolling_displacement<double>(math::as_sphere(*foot), Eigen::Quaterniond::Identity(),
                    Eigen::Quaterniond(Eigen::AngleAxisd(0.1, Eigen::Vector3d::UnitX())), Eigen::Vector3d::UnitZ());
    EXPECT_TRUE(displacement.isApprox(Eigen::Vector3d(0.0, -0.0033, 0.0)));
}

TEST(rolling_shape, type_and_as_shape_recover_the_shape) {
    const math::Sphere<double> sphere(0.033);
    const math::Torus<double> torus(Eigen::Vector3d(0.1, 1.0, -0.2), 0.1, 0.02);
    const math::RollingShape<double>& sphere_shape = sphere;
    const math::RollingShape<double>& torus_shape = torus;
    EXPECT_TRUE(sphere_shape.type() == math::RollingShapeType::SPHERE);
    EXPECT_TRUE(torus_shape.type() == math::RollingShapeType::TORUS);
    EXPECT_EQ(&math::as_sphere(sphere_shape), &sphere);
    EXPECT_EQ(&math::as_torus(torus_shape), &torus);
    EXPECT_THROW(math::as_torus(sphere_shape), std::bad_cast);
    EXPECT_THROW(math::as_sphere(torus_shape), std::bad_cast);
}

TEST(rolling_shape, is_rolling_shape) {
    EXPECT_TRUE(math::IsRollingShape<math::Sphere<double>>);
    EXPECT_TRUE(math::IsRollingShape<math::Torus<float>>);
    // A type with a scalar that is not derived from RollingShape, and a type without one
    EXPECT_FALSE(math::IsRollingShape<Eigen::Vector3d>);
    EXPECT_FALSE(math::IsRollingShape<double>);
}

TEST(rolling_shape, rolling_displacement_is_independent_of_the_signs_of_the_orientations) {
    const math::Torus<double> wheel(Eigen::Vector3d(0.1, 1.0, -0.2), 0.1, 0.02);
    const Eigen::Vector3d normal = Eigen::Vector3d(0.1, 0.2, 1.0).normalized();
    const Eigen::Quaterniond start_orientation(Eigen::AngleAxisd(0.7, Eigen::Vector3d(1.0, 2.0, 3.0).normalized()));
    const Eigen::Quaterniond end_orientation =
            Eigen::Quaterniond(Eigen::AngleAxisd(0.3, Eigen::Vector3d(-0.5, 1.0, 0.2).normalized())) *
            start_orientation;
    const Eigen::Vector3d displacement =
            math::rolling_displacement<double>(wheel, start_orientation, end_orientation, normal);
    EXPECT_TRUE(math::rolling_displacement<double>(wheel, start_orientation,
            Eigen::Quaterniond(-end_orientation.coeffs()), normal)
                    .isApprox(displacement));
    EXPECT_TRUE(math::rolling_displacement<double>(wheel, Eigen::Quaterniond(-start_orientation.coeffs()),
            end_orientation, normal)
                    .isApprox(displacement));
}

TEST(rolling_shape, rolling_displacement_through_orientations) {
    const math::Sphere<double> sphere(0.033);
    const Eigen::Vector3d normal = Eigen::Vector3d::UnitZ();
    const std::vector<Eigen::Quaterniond> orientations{Eigen::Quaterniond::Identity(),
            Eigen::Quaterniond(Eigen::AngleAxisd(0.1, Eigen::Vector3d::UnitX())),
            Eigen::Quaterniond(Eigen::AngleAxisd(0.2, Eigen::Vector3d(1.0, 1.0, 0.0).normalized()))};
    EXPECT_TRUE(math::rolling_displacement<double>(sphere, orientations, normal)
                    .isApprox(math::rolling_displacement<double>(sphere, orientations[0], orientations[1], normal) +
                              math::rolling_displacement<double>(sphere, orientations[1], orientations[2], normal)));
    EXPECT_TRUE(math::rolling_displacement<double>(sphere, std::vector<Eigen::Quaterniond>{orientations[1]}, normal)
                    .isZero());
    EXPECT_THROW(math::rolling_displacement<double>(sphere, std::vector<Eigen::Quaterniond>{}, normal),
            std::invalid_argument);
}

TEST(sphere, documentation_example) {
    const math::Sphere<double> foot(0.033);
    const Eigen::Vector3d contact_to_centre =
            foot.contact_to_centre<double>(Eigen::Quaterniond::Identity(), Eigen::Vector3d::UnitZ());
    EXPECT_TRUE(contact_to_centre.isApprox(Eigen::Vector3d(0.0, 0.0, 0.033)));
}

TEST(sphere, contact_to_centre_is_the_radius_along_the_normal) {
    const math::Sphere<double> sphere(0.033);
    const Eigen::Vector3d normal = Eigen::Vector3d(0.1, -0.2, 1.0).normalized();
    for (int i = 0; i < 10; ++i) {
        EXPECT_TRUE(
                sphere.contact_to_centre<double>(Eigen::Quaterniond::UnitRandom(), normal).isApprox(0.033 * normal));
    }
}

TEST(sphere, rotating_about_a_fixed_axis) {
    const double radius = 0.033, angle = 0.8;
    const math::Sphere<double> sphere(radius);
    const Eigen::Vector3d normal = Eigen::Vector3d(0.2, -0.1, 1.0).normalized();
    const Eigen::Vector3d axis = Eigen::Vector3d(0.3, 1.0, 0.4).normalized();
    const Eigen::Quaterniond start_orientation = Eigen::Quaterniond::UnitRandom();
    EXPECT_TRUE(math::rolling_displacement<double>(sphere, start_orientation,
            Eigen::Quaterniond(Eigen::AngleAxisd(angle, axis)) * start_orientation, normal)
                    .isApprox(radius * angle * axis.cross(normal)));
}

TEST(sphere, spinning_about_the_normal_does_not_move) {
    const math::Sphere<double> sphere(0.033);
    const Eigen::Vector3d normal = Eigen::Vector3d(0.2, -0.1, 1.0).normalized();
    const Eigen::Quaterniond start_orientation = Eigen::Quaterniond::UnitRandom();
    EXPECT_LT(math::rolling_displacement<double>(sphere, start_orientation,
                      Eigen::Quaterniond(Eigen::AngleAxisd(0.6, normal)) * start_orientation, normal)
                      .norm(),
            1.0e-15);
}

TEST(sphere, increments_follow_the_rolling_path) {
    check_increments_follow_the_rolling_path(math::Sphere<double>(0.033));
}

TEST(sphere, negative_radius_throws) {
    EXPECT_THROW(math::Sphere<double>(-0.01), std::invalid_argument);
}

TEST(torus, documentation_example) {
    const math::Torus<double> wheel(Eigen::Vector3d::UnitY(), 0.1, 0.02);
    const Eigen::Vector3d displacement = math::rolling_displacement<double>(wheel, Eigen::Quaterniond::Identity(),
            Eigen::Quaterniond(Eigen::AngleAxisd(0.5, Eigen::Vector3d::UnitY())), Eigen::Vector3d::UnitZ());
    EXPECT_TRUE(displacement.isApprox(Eigen::Vector3d(0.06, 0.0, 0.0)));
}

TEST(torus, axis_is_normalised) {
    EXPECT_TRUE(
            math::Torus<double>(Eigen::Vector3d(0.0, 2.0, 0.0), 0.1, 0.0).axis().isApprox(Eigen::Vector3d::UnitY()));
}

TEST(torus, upright_contact_to_centre_is_both_radii_along_the_normal) {
    const math::Torus<double> wheel(Eigen::Vector3d::UnitY(), 0.1, 0.02);
    // Any heading and any spin about the axis keep the axis horizontal
    const Eigen::Quaterniond orientation = Eigen::Quaterniond(Eigen::AngleAxisd(0.7, Eigen::Vector3d::UnitZ())) *
                                           Eigen::Quaterniond(Eigen::AngleAxisd(1.3, Eigen::Vector3d::UnitY()));
    EXPECT_TRUE(wheel.contact_to_centre<double>(orientation, Eigen::Vector3d::UnitZ())
                    .isApprox(0.12 * Eigen::Vector3d::UnitZ()));
}

TEST(torus, cambered_contact_to_centre) {
    const double major_radius = 0.1, minor_radius = 0.02, camber = 0.2;
    const math::Torus<double> wheel(Eigen::Vector3d::UnitY(), major_radius, minor_radius);
    // Camber about x tilts the axis to (0, cos, sin), and the major circle's lowest point is below the centre along
    // (0, -sin, cos)
    EXPECT_TRUE(
            wheel.contact_to_centre<double>(Eigen::Quaterniond(Eigen::AngleAxisd(camber, Eigen::Vector3d::UnitX())),
                         Eigen::Vector3d::UnitZ())
                    .isApprox(major_radius * Eigen::Vector3d(0.0, -std::sin(camber), std::cos(camber)) +
                              minor_radius * Eigen::Vector3d::UnitZ()));
}

TEST(torus, cambered_rolling_about_its_axis) {
    const double major_radius = 0.1, minor_radius = 0.02, camber = 0.2, heading = 0.4, angle = 0.5;
    const math::Torus<double> wheel(Eigen::Vector3d::UnitY(), major_radius, minor_radius);
    const Eigen::Quaterniond heading_rotation(Eigen::AngleAxisd(heading, Eigen::Vector3d::UnitZ()));
    const Eigen::Quaterniond start_orientation =
            heading_rotation * Eigen::Quaterniond(Eigen::AngleAxisd(camber, Eigen::Vector3d::UnitX()));
    const Eigen::Vector3d axis = start_orientation * Eigen::Vector3d::UnitY();
    // Cambered direction from the major circle's lowest point to the centre turned by the heading, plus the minor
    // radius along the normal
    const Eigen::Vector3d contact_to_centre =
            major_radius * (heading_rotation * Eigen::Vector3d(0.0, -std::sin(camber), std::cos(camber))) +
            minor_radius * Eigen::Vector3d::UnitZ();
    EXPECT_TRUE(math::rolling_displacement<double>(wheel, start_orientation,
            Eigen::Quaterniond(Eigen::AngleAxisd(angle, axis)) * start_orientation, Eigen::Vector3d::UnitZ())
                    .isApprox(angle * axis.cross(contact_to_centre)));
}

TEST(torus, increments_follow_the_rolling_path) {
    check_increments_follow_the_rolling_path(math::Torus<double>(Eigen::Vector3d(0.1, 1.0, -0.2), 0.1, 0.02));
}

TEST(torus, lying_flat_throws) {
    const math::Torus<double> wheel(Eigen::Vector3d::UnitZ(), 0.1, 0.02);
    EXPECT_THROW(wheel.contact_to_centre<double>(Eigen::Quaterniond::Identity(), Eigen::Vector3d::UnitZ()),
            std::invalid_argument);
}

TEST(torus, invalid_parameters_throw) {
    EXPECT_THROW(math::Torus<double>(Eigen::Vector3d::Zero(), 0.1, 0.02), std::invalid_argument);
    EXPECT_THROW(math::Torus<double>(Eigen::Vector3d::UnitY(), 0.0, 0.02), std::invalid_argument);
    EXPECT_THROW(math::Torus<double>(Eigen::Vector3d::UnitY(), 0.1, -0.01), std::invalid_argument);
}
