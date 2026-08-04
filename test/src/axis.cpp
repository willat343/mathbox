#include "mathbox/axis.hpp"

#include <gtest/gtest.h>

#include <numbers>

TEST(AxisType, values) {
    EXPECT_EQ(math::AxisType::count, 3u);
    EXPECT_EQ(math::AxisType::X, 0u);
    EXPECT_EQ(math::AxisType::Y, 1u);
    EXPECT_EQ(math::AxisType::Z, 2u);
}

TEST(SignedAxis, positive_for_and_negative_for) {
    EXPECT_EQ(math::SignedAxis::positive_for(math::AxisType::X), math::SignedAxisType::POSITIVE_X);
    EXPECT_EQ(math::SignedAxis::negative_for(math::AxisType::X), math::SignedAxisType::NEGATIVE_X);
    EXPECT_EQ(math::SignedAxis::positive_for(math::AxisType::Y), math::SignedAxisType::POSITIVE_Y);
    EXPECT_EQ(math::SignedAxis::negative_for(math::AxisType::Y), math::SignedAxisType::NEGATIVE_Y);
    EXPECT_EQ(math::SignedAxis::positive_for(math::AxisType::Z), math::SignedAxisType::POSITIVE_Z);
    EXPECT_EQ(math::SignedAxis::negative_for(math::AxisType::Z), math::SignedAxisType::NEGATIVE_Z);
}

TEST(SignedAxis, axis) {
    EXPECT_EQ(math::SignedAxis(math::SignedAxisType::POSITIVE_X).axis(), math::AxisType::X);
    EXPECT_EQ(math::SignedAxis(math::SignedAxisType::NEGATIVE_X).axis(), math::AxisType::X);
    EXPECT_EQ(math::SignedAxis(math::SignedAxisType::POSITIVE_Y).axis(), math::AxisType::Y);
    EXPECT_EQ(math::SignedAxis(math::SignedAxisType::NEGATIVE_Y).axis(), math::AxisType::Y);
    EXPECT_EQ(math::SignedAxis(math::SignedAxisType::POSITIVE_Z).axis(), math::AxisType::Z);
    EXPECT_EQ(math::SignedAxis(math::SignedAxisType::NEGATIVE_Z).axis(), math::AxisType::Z);
}

TEST(SignedAxis, is_positive_and_is_negative) {
    const math::SignedAxis positive_x(math::SignedAxisType::POSITIVE_X);
    const math::SignedAxis negative_x(math::SignedAxisType::NEGATIVE_X);
    EXPECT_TRUE(positive_x.is_positive());
    EXPECT_FALSE(positive_x.is_negative());
    EXPECT_FALSE(negative_x.is_positive());
    EXPECT_TRUE(negative_x.is_negative());
}

TEST(SignedAxis, positive_and_negative) {
    const math::SignedAxis positive_y(math::SignedAxisType::POSITIVE_Y);
    const math::SignedAxis negative_y(math::SignedAxisType::NEGATIVE_Y);
    EXPECT_EQ(positive_y.positive(), positive_y);
    EXPECT_EQ(positive_y.negative(), negative_y);
    EXPECT_EQ(negative_y.positive(), positive_y);
    EXPECT_EQ(negative_y.negative(), negative_y);
}

TEST(SignedAxis, sign) {
    EXPECT_DOUBLE_EQ((math::SignedAxis(math::SignedAxisType::POSITIVE_Z).sign<double>()), 1.0);
    EXPECT_DOUBLE_EQ((math::SignedAxis(math::SignedAxisType::NEGATIVE_Z).sign<double>()), -1.0);
}

TEST(SignedAxis, direction) {
    EXPECT_TRUE((
            math::SignedAxis(math::SignedAxisType::POSITIVE_X).direction<double>().isApprox(Eigen::Vector3d::UnitX())));
    EXPECT_TRUE((math::SignedAxis(math::SignedAxisType::NEGATIVE_X)
                    .direction<double>()
                    .isApprox(-Eigen::Vector3d::UnitX())));
    EXPECT_TRUE((
            math::SignedAxis(math::SignedAxisType::POSITIVE_Y).direction<double>().isApprox(Eigen::Vector3d::UnitY())));
    EXPECT_TRUE((math::SignedAxis(math::SignedAxisType::NEGATIVE_Y)
                    .direction<double>()
                    .isApprox(-Eigen::Vector3d::UnitY())));
    EXPECT_TRUE((
            math::SignedAxis(math::SignedAxisType::POSITIVE_Z).direction<double>().isApprox(Eigen::Vector3d::UnitZ())));
    EXPECT_TRUE((math::SignedAxis(math::SignedAxisType::NEGATIVE_Z)
                    .direction<double>()
                    .isApprox(-Eigen::Vector3d::UnitZ())));
}

TEST(SignedAxis, for_vector) {
    EXPECT_EQ(math::SignedAxis::for_vector(Eigen::Vector3d(9.81, 0.1, -0.2)), math::SignedAxisType::POSITIVE_X);
    EXPECT_EQ(math::SignedAxis::for_vector(Eigen::Vector3d(-9.81, 0.1, -0.2)), math::SignedAxisType::NEGATIVE_X);
    EXPECT_EQ(math::SignedAxis::for_vector(Eigen::Vector3d(0.1, 9.81, -0.2)), math::SignedAxisType::POSITIVE_Y);
    EXPECT_EQ(math::SignedAxis::for_vector(Eigen::Vector3d(0.1, -9.81, -0.2)), math::SignedAxisType::NEGATIVE_Y);
    EXPECT_EQ(math::SignedAxis::for_vector(Eigen::Vector3d(0.1, -0.2, 9.81)), math::SignedAxisType::POSITIVE_Z);
    EXPECT_EQ(math::SignedAxis::for_vector(Eigen::Vector3d(0.1, -0.2, -9.81)), math::SignedAxisType::NEGATIVE_Z);
}

TEST(SignedAxis, angle_to) {
    const math::SignedAxis positive_z(math::SignedAxisType::POSITIVE_Z);
    EXPECT_NEAR(positive_z.angle_to(Eigen::Vector3d::UnitZ()), 0.0, 1.0e-12);
    EXPECT_NEAR(positive_z.angle_to(Eigen::Vector3d(1.0, 0.0, 1.0)), std::numbers::pi / 4.0, 1.0e-12);
    EXPECT_NEAR(positive_z.angle_to(-Eigen::Vector3d::UnitZ()), std::numbers::pi, 1.0e-12);
}
