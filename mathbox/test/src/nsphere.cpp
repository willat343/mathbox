#include "mathbox/nsphere.hpp"

#include <gtest/gtest.h>

TEST(nsphere, coordinates_1D_0) {
    Eigen::Matrix<double, 1, 1> cartesian_coordinates(0.0);
    EXPECT_TRUE(cartesian_coordinates.isApprox(
            math::spherical_to_cartesian(math::cartesian_to_spherical(cartesian_coordinates))));
}

TEST(nsphere, coordinates_1D_1) {
    Eigen::Matrix<double, 1, 1> cartesian_coordinates(1.0);
    EXPECT_TRUE(cartesian_coordinates.isApprox(
            math::spherical_to_cartesian(math::cartesian_to_spherical(cartesian_coordinates))));
}

TEST(nsphere, coordinates_1D_2) {
    Eigen::Matrix<double, 1, 1> cartesian_coordinates(1.7);
    EXPECT_TRUE(cartesian_coordinates.isApprox(
            math::spherical_to_cartesian(math::cartesian_to_spherical(cartesian_coordinates))));
}

TEST(nsphere, coordinates_1D_3) {
    Eigen::Matrix<double, 1, 1> cartesian_coordinates(-2.1);
    EXPECT_TRUE(cartesian_coordinates.isApprox(
            math::spherical_to_cartesian(math::cartesian_to_spherical(cartesian_coordinates))));
}

TEST(nsphere, coordinates_2D_0) {
    Eigen::Vector2d cartesian_coordinates(0.0, 0.0);
    EXPECT_TRUE(cartesian_coordinates.isApprox(
            math::spherical_to_cartesian(math::cartesian_to_spherical(cartesian_coordinates))));
}

TEST(nsphere, coordinates_2D_1) {
    Eigen::Vector2d cartesian_coordinates(1.0, 0.0);
    EXPECT_TRUE(cartesian_coordinates.isApprox(
            math::spherical_to_cartesian(math::cartesian_to_spherical(cartesian_coordinates))));
    Eigen::Vector2d unit_cartesian_coordinates = cartesian_coordinates.normalized();
    EXPECT_TRUE(unit_cartesian_coordinates.isApprox(math::spherical_angles_to_unit_cartesian<2>(
            math::cartesian_to_spherical_angles(unit_cartesian_coordinates))));
    EXPECT_TRUE(math::cartesian_to_spherical(cartesian_coordinates)
                        .tail<1>()
                        .isApprox(math::cartesian_to_spherical_angles(unit_cartesian_coordinates)));
}

TEST(nsphere, coordinates_2D_2) {
    Eigen::Vector2d cartesian_coordinates(1.7, 0.6);
    EXPECT_TRUE(cartesian_coordinates.isApprox(
            math::spherical_to_cartesian(math::cartesian_to_spherical(cartesian_coordinates))));
    Eigen::Vector2d unit_cartesian_coordinates = cartesian_coordinates.normalized();
    EXPECT_TRUE(unit_cartesian_coordinates.isApprox(math::spherical_angles_to_unit_cartesian<2>(
            math::cartesian_to_spherical_angles(unit_cartesian_coordinates))));
    EXPECT_TRUE(math::cartesian_to_spherical(cartesian_coordinates)
                        .tail<1>()
                        .isApprox(math::cartesian_to_spherical_angles(unit_cartesian_coordinates)));
}

TEST(nsphere, coordinates_2D_3) {
    Eigen::Vector2d cartesian_coordinates(-2.1, -3.5);
    EXPECT_TRUE(cartesian_coordinates.isApprox(
            math::spherical_to_cartesian(math::cartesian_to_spherical(cartesian_coordinates))));
    Eigen::Vector2d unit_cartesian_coordinates = cartesian_coordinates.normalized();
    EXPECT_TRUE(unit_cartesian_coordinates.isApprox(math::spherical_angles_to_unit_cartesian<2>(
            math::cartesian_to_spherical_angles(unit_cartesian_coordinates))));
    EXPECT_TRUE(math::cartesian_to_spherical(cartesian_coordinates)
                        .tail<1>()
                        .isApprox(math::cartesian_to_spherical_angles(unit_cartesian_coordinates)));
}

TEST(nsphere, coordinates_3D_0) {
    Eigen::Vector3d cartesian_coordinates(0.0, 0.0, 0.0);
    EXPECT_TRUE(cartesian_coordinates.isApprox(
            math::spherical_to_cartesian(math::cartesian_to_spherical(cartesian_coordinates))));
}

TEST(nsphere, coordinates_3D_1) {
    Eigen::Vector3d cartesian_coordinates(1.0, 0.0, 0.0);
    EXPECT_TRUE(cartesian_coordinates.isApprox(
            math::spherical_to_cartesian(math::cartesian_to_spherical(cartesian_coordinates))));
    Eigen::Vector3d unit_cartesian_coordinates = cartesian_coordinates.normalized();
    EXPECT_TRUE(unit_cartesian_coordinates.isApprox(math::spherical_angles_to_unit_cartesian<3>(
            math::cartesian_to_spherical_angles(unit_cartesian_coordinates))));
    EXPECT_TRUE(math::cartesian_to_spherical(cartesian_coordinates)
                        .tail<2>()
                        .isApprox(math::cartesian_to_spherical_angles(unit_cartesian_coordinates)));
}

TEST(nsphere, coordinates_3D_2) {
    Eigen::Vector3d cartesian_coordinates(1.7, 0.6, -0.7);
    EXPECT_TRUE(cartesian_coordinates.isApprox(
            math::spherical_to_cartesian(math::cartesian_to_spherical(cartesian_coordinates))));
    Eigen::Vector3d unit_cartesian_coordinates = cartesian_coordinates.normalized();
    EXPECT_TRUE(unit_cartesian_coordinates.isApprox(math::spherical_angles_to_unit_cartesian<3>(
            math::cartesian_to_spherical_angles(unit_cartesian_coordinates))));
    EXPECT_TRUE(math::cartesian_to_spherical(cartesian_coordinates)
                        .tail<2>()
                        .isApprox(math::cartesian_to_spherical_angles(unit_cartesian_coordinates)));
}

TEST(nsphere, coordinates_3D_3) {
    Eigen::Vector3d cartesian_coordinates(-2.1, -3.5, -6.6);
    EXPECT_TRUE(cartesian_coordinates.isApprox(
            math::spherical_to_cartesian(math::cartesian_to_spherical(cartesian_coordinates))));
    Eigen::Vector3d unit_cartesian_coordinates = cartesian_coordinates.normalized();
    EXPECT_TRUE(unit_cartesian_coordinates.isApprox(math::spherical_angles_to_unit_cartesian<3>(
            math::cartesian_to_spherical_angles(unit_cartesian_coordinates))));
    EXPECT_TRUE(math::cartesian_to_spherical(cartesian_coordinates)
                        .tail<2>()
                        .isApprox(math::cartesian_to_spherical_angles(unit_cartesian_coordinates)));
}
