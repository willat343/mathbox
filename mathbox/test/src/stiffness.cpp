#include "mathbox/stiffness.hpp"

#include <gtest/gtest.h>

#include <Eigen/Core>

TEST(stiffness, I) {
    EXPECT_EQ((math::stiffness_from_covariance<3, double>(Eigen::Matrix3d::Identity())), Eigen::Matrix3d::Identity());
    EXPECT_EQ((math::stiffness_from_information<3, double>(Eigen::Matrix3d::Identity())), Eigen::Matrix3d::Identity());
}

TEST(stiffness, ascending) {
    const Eigen::Matrix3d covariance = (Eigen::Vector3d() << 1.0, 2.0, 3.0).finished().asDiagonal();
    const Eigen::Matrix3d expected_stiffness = covariance.inverse().cwiseSqrt();
    EXPECT_TRUE(math::stiffness_from_covariance(covariance).isApprox(expected_stiffness));
}

TEST(stiffness, descending) {
    const Eigen::Matrix3d covariance = (Eigen::Vector3d() << 3.0, 2.0, 1.0).finished().asDiagonal();
    const Eigen::Matrix3d expected_stiffness = covariance.inverse().cwiseSqrt();
    EXPECT_TRUE(math::stiffness_from_covariance(covariance).isApprox(expected_stiffness));
}

TEST(stiffness, mixed_1) {
    const Eigen::Matrix3d covariance = (Eigen::Vector3d() << 1.0, 3.0, 2.0).finished().asDiagonal();
    const Eigen::Matrix3d expected_stiffness = covariance.inverse().cwiseSqrt();
    EXPECT_TRUE(math::stiffness_from_covariance(covariance).isApprox(expected_stiffness));
}

TEST(stiffness, mixed_2) {
    const Eigen::Matrix3d covariance = (Eigen::Vector3d() << 2.0, 1.0, 3.0).finished().asDiagonal();
    const Eigen::Matrix3d expected_stiffness = covariance.inverse().cwiseSqrt();
    EXPECT_TRUE(math::stiffness_from_covariance(covariance).isApprox(expected_stiffness));
}

TEST(stiffness, mixed_3) {
    const Eigen::Matrix3d covariance = (Eigen::Vector3d() << 3.0, 1.0, 2.0).finished().asDiagonal();
    const Eigen::Matrix3d expected_stiffness = covariance.inverse().cwiseSqrt();
    EXPECT_TRUE(math::stiffness_from_covariance(covariance).isApprox(expected_stiffness));
}

TEST(stiffness, mixed_4) {
    const Eigen::Matrix3d covariance = (Eigen::Vector3d() << 2.0, 3.0, 1.0).finished().asDiagonal();
    const Eigen::Matrix3d expected_stiffness = covariance.inverse().cwiseSqrt();
    EXPECT_TRUE(math::stiffness_from_covariance(covariance).isApprox(expected_stiffness));
}

TEST(stiffness, symmetric_positive_definite) {
    const Eigen::Matrix3d covariance = (Eigen::Matrix3d() << 2.5, 0.4, -0.9, 0.4, 0.8, 0.9, -0.9, 0.9, 6.3).finished();
    const Eigen::Matrix3d stiffness = math::stiffness_from_covariance(covariance);
    EXPECT_TRUE(covariance.inverse().isApprox(stiffness.transpose() * stiffness));
}
