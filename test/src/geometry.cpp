#include "mathbox/geometry.hpp"

#include <gtest/gtest.h>

#include "mathbox/matrix_operations.hpp"

void check_transform_adjoint_blocks(const Eigen::Matrix<double, 6, 6>& transform_adjoint,
        const Eigen::Isometry3d& transform) {
    const Eigen::Matrix3d R = transform.rotation();
    const Eigen::Vector3d t = transform.translation();
    const Eigen::Matrix3d t_SS = math::skew_symmetric_cross(t);
    const Eigen::Matrix3d t_SS_times_R = t_SS * R;
    const Eigen::Matrix3d top_left = transform_adjoint.block<3, 3>(0, 0);
    const Eigen::Matrix3d top_right = transform_adjoint.block<3, 3>(0, 3);
    const Eigen::Matrix3d bottom_left = transform_adjoint.block<3, 3>(3, 0);
    const Eigen::Matrix3d bottom_right = transform_adjoint.block<3, 3>(3, 3);
    EXPECT_TRUE(top_left.isApprox(R));
    EXPECT_TRUE(top_right.isApprox(Eigen::Matrix3d::Zero()));
    EXPECT_TRUE(bottom_left.isApprox(t_SS_times_R));
    EXPECT_TRUE(bottom_right.isApprox(R));
}

TEST(change_relative_transform_frame, identity_with_identity) {
    const Eigen::Isometry3d relative_transform_A = Eigen::Isometry3d::Identity();
    const Eigen::Isometry3d transform_B_A =
            Eigen::Translation<double, 3>{1.0, 2.0, 3.0} * Eigen::Quaterniond{0.8, 0.1, 0.05, 0.2}.normalized();
    const Eigen::Isometry3d relative_transform_B =
            math::change_relative_transform_frame(relative_transform_A, transform_B_A);
    const Eigen::Isometry3d expected = Eigen::Isometry3d::Identity();
    EXPECT_TRUE(relative_transform_B.isApprox(expected));
}

TEST(change_relative_transform_frame, transform_with_identity) {
    const Eigen::Isometry3d relative_transform_A =
            Eigen::Translation<double, 3>{-5.0, 6.0, -10.0} * Eigen::Quaterniond{0.4, 0.2, 0.7, 0.15}.normalized();
    const Eigen::Isometry3d transform_B_A = Eigen::Isometry3d::Identity();
    const Eigen::Isometry3d relative_transform_B =
            math::change_relative_transform_frame(relative_transform_A, transform_B_A);
    const Eigen::Isometry3d expected = relative_transform_A;
    EXPECT_TRUE(relative_transform_B.isApprox(expected));
}

TEST(change_relative_transform_frame, translation_with_translation) {
    const Eigen::Isometry3d relative_transform_A{Eigen::Translation<double, 3>{-5.0, 6.0, -10.0}};
    const Eigen::Isometry3d transform_B_A{Eigen::Translation<double, 3>{1.0, 2.0, 3.0}};
    const Eigen::Isometry3d relative_transform_B =
            math::change_relative_transform_frame(relative_transform_A, transform_B_A);
    const Eigen::Isometry3d expected = relative_transform_A;
    EXPECT_TRUE(relative_transform_B.isApprox(expected));
}

TEST(change_relative_transform_frame, rotation_with_rotation) {
    const Eigen::Isometry3d relative_transform_A{Eigen::Quaterniond{0.4, 0.2, 0.7, 0.15}.normalized()};
    const Eigen::Isometry3d transform_B_A{Eigen::Quaterniond{0.8, 0.1, 0.05, 0.2}.normalized()};
    const Eigen::Isometry3d relative_transform_B =
            math::change_relative_transform_frame(relative_transform_A, transform_B_A);
    EXPECT_NEAR(Eigen::AngleAxisd(relative_transform_B.rotation()).angle(),
            Eigen::AngleAxisd(relative_transform_A.rotation()).angle(), 1.0e-12);
}

TEST(change_relative_transform_frame, check_against_alternate_method) {
    const Eigen::Quaterniond q = Eigen::Quaterniond{0.7, 0.123, -0.2, -0.62}.normalized();
    const Eigen::Translation3d t{5.6, 9.2, -5.5};
    const Eigen::Isometry3d T_rel = t * q;
    const Eigen::Quaterniond q_r = Eigen::Quaterniond{0.1, 0.55, 0.8, -0.622}.normalized();
    const Eigen::Translation3d t_r{-15.1, 0.55, 2.8};
    const Eigen::Isometry3d T_rigid = t_r * q_r;
    const Eigen::Isometry3d T_method = math::change_relative_transform_frame(T_rel, T_rigid);
    const Eigen::Isometry3d T_alt = (T_rigid * (T_rigid * T_rel).inverse()).inverse();
    EXPECT_TRUE(T_method.isApprox(T_alt));
}

TEST(change_tf_covariance_frame, identity) {
    Eigen::Matrix<double, 6, 6> covariance;
    covariance << 0.1, 0.02, 0.03, 0.04, 0.05, 0.06, 0.02, 0.2, 0.07, 0.08, 0.09, 0.10, 0.03, 0.07, 0.3, 0.11, 0.12,
            0.13, 0.04, 0.08, 0.11, 0.4, 0.14, 0.15, 0.05, 0.09, 0.12, 0.14, 0.5, 0.16, 0.06, 0.10, 0.13, 0.15, 0.16,
            0.6;
    const Eigen::Isometry3d transform = Eigen::Isometry3d::Identity();
    EXPECT_TRUE(covariance.isApprox(math::change_tf_covariance_frame(covariance, transform)));
}

TEST(compute_constant_rates, compute_constant_rates_0) {
    const Eigen::Isometry3d I = Eigen::Isometry3d::Identity();
    const Eigen::Vector3d rotation_axis = Eigen::Vector3d(0.8, 0.1, 0.05).normalized();
    const double rotation_angle_rate{0.492};
    const Eigen::Vector3d linear_velocity{1.0, 2.0, 3.0};
    const double dt = 1.0;
    const Eigen::Vector3d angular_velocity = rotation_angle_rate * dt * rotation_axis;
    const Eigen::Isometry3d pose = Eigen::Translation<double, 3>{linear_velocity * dt} *
                                   Eigen::AngleAxisd{rotation_angle_rate * dt, rotation_axis};
    Eigen::Matrix<double, 6, 1> rates = math::compute_constant_rates(I, pose, dt);
    const Eigen::Vector3d angular_velocity_out = rates.block<3, 1>(0, 0);
    const Eigen::Vector3d linear_velocity_out = rates.block<3, 1>(3, 0);
    EXPECT_TRUE(angular_velocity_out.isApprox(angular_velocity));
    EXPECT_TRUE(linear_velocity_out.isApprox(linear_velocity));
}

TEST(compute_constant_rates, compute_constant_rates_1) {
    const Eigen::Isometry3d I = Eigen::Isometry3d::Identity();
    const Eigen::Vector3d rotation_axis = Eigen::Vector3d(0.8, 0.1, 0.05).normalized();
    const double rotation_angle_rate{0.492};
    const Eigen::Vector3d linear_velocity{1.0, 2.0, 3.0};
    const double dt = M_PI / rotation_angle_rate;  // Too large a dt will invalidate test because the rotation will wrap
    const Eigen::Vector3d angular_velocity = rotation_angle_rate * rotation_axis;
    const Eigen::Isometry3d pose = Eigen::Translation<double, 3>{linear_velocity * dt} *
                                   Eigen::AngleAxisd{rotation_angle_rate * dt, rotation_axis};
    Eigen::Matrix<double, 6, 1> rates = math::compute_constant_rates(I, pose, dt);
    const Eigen::Vector3d angular_velocity_out = rates.block<3, 1>(0, 0);
    const Eigen::Vector3d linear_velocity_out = rates.block<3, 1>(3, 0);
    EXPECT_TRUE(angular_velocity_out.isApprox(angular_velocity));
    EXPECT_TRUE(linear_velocity_out.isApprox(linear_velocity));
}

TEST(compute_constant_rates, compute_constant_rates_2) {
    const Eigen::Vector3d start_position{0.5, -5.0, 0.0};
    const Eigen::Isometry3d start_pose = Eigen::Translation<double, 3>{start_position} * Eigen::Quaterniond::Identity();
    const Eigen::Vector3d end_position{2.0, 3.0, 4.0};
    const Eigen::Isometry3d end_pose = Eigen::Translation<double, 3>{2.0, 3.0, 4.0} * Eigen::Quaterniond::Identity();
    const double dt{5.0};
    const Eigen::Matrix<double, 6, 1> rates = math::compute_constant_rates(start_pose, end_pose, dt);
    Eigen::Matrix<double, 6, 1> rates_reference;
    rates_reference << Eigen::Vector3d::Zero(), (end_position - start_position) / dt;
    EXPECT_TRUE(rates.isApprox(rates_reference));
}

TEST(glerp, identity_identity) {
    const Eigen::Isometry3d T_0 = Eigen::Isometry3d::Identity();
    const Eigen::Isometry3d T_1 = Eigen::Isometry3d::Identity();
    EXPECT_TRUE(math::glerp(T_0, T_1, 0.0).isApprox(T_0));
    EXPECT_TRUE(math::glerp(T_0, T_1, 1.0).isApprox(T_1));
}

TEST(glerp, lhs_identity) {
    const Eigen::Isometry3d T_0 =
            Eigen::Translation<double, 3>{1.0, 2.0, 3.0} * Eigen::Quaterniond{0.8, 0.1, 0.05, 0.2}.normalized();
    const Eigen::Isometry3d T_1 = Eigen::Isometry3d::Identity();
    EXPECT_TRUE(math::glerp(T_0, T_1, 0.0).isApprox(T_0));
    EXPECT_TRUE(math::glerp(T_0, T_1, 1.0).isApprox(T_1));
}

TEST(glerp, rhs_identity) {
    const Eigen::Isometry3d T_0 = Eigen::Isometry3d::Identity();
    const Eigen::Isometry3d T_1 =
            Eigen::Translation<double, 3>{1.0, 2.0, 3.0} * Eigen::Quaterniond{0.8, 0.1, 0.05, 0.2}.normalized();
    EXPECT_TRUE(math::glerp(T_0, T_1, 0.0).isApprox(T_0));
    EXPECT_TRUE(math::glerp(T_0, T_1, 1.0).isApprox(T_1));
}

TEST(relative_transform, transform_0) {
    const Eigen::Isometry3d I = Eigen::Isometry3d::Identity();
    const Eigen::Isometry3d pose =
            Eigen::Translation<double, 3>{1.0, 2.0, 3.0} * Eigen::Quaterniond{0.8, 0.1, 0.05, 0.2}.normalized();
    const Eigen::Isometry3d transform = math::relative_transform(I, pose);
    EXPECT_TRUE(transform.isApprox(pose));
}

TEST(relative_transform, transform_0_inv) {
    const Eigen::Isometry3d I = Eigen::Isometry3d::Identity();
    const Eigen::Vector3d translation{1.0, 2.0, 3.0};
    const Eigen::Isometry3d pose =
            Eigen::Translation<double, 3>{translation} * Eigen::Quaterniond{0.8, 0.1, 0.05, 0.2}.normalized();
    const Eigen::Isometry3d transform = math::relative_transform(pose, I);
    EXPECT_TRUE(transform.isApprox(pose.inverse()));
}

TEST(relative_transform_change_frame, identity_with_transform) {
    const Eigen::Isometry3d relative_transform_A = Eigen::Isometry3d::Identity();
    const Eigen::Isometry3d transform_B_A = Eigen::Isometry3d::Identity();
    const Eigen::Isometry3d relative_transform_B =
            math::change_relative_transform_frame(relative_transform_A, transform_B_A);
    const Eigen::Isometry3d expected = Eigen::Isometry3d::Identity();
    EXPECT_TRUE(relative_transform_B.isApprox(expected));
}

TEST(rotate_point_covariance, identity_rotation) {
    Eigen::Matrix<double, 3, 3> covariance;
    covariance << 0.1, 0.002, 0.003, 0.002, 0.2, 0.006, 0.003, 0.006, 0.3;
    Eigen::Matrix3d rotation = Eigen::Matrix3d::Identity();
    EXPECT_TRUE(covariance.isApprox(math::rotate_point_covariance(covariance, rotation)));
}

TEST(rotate_point_covariance, yaw_90_x_covariance) {
    Eigen::Matrix<double, 3, 3> covariance;
    covariance << 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    Eigen::AngleAxisd rotation = Eigen::AngleAxisd(M_PI / 2.0, Eigen::Vector3d::UnitZ());
    Eigen::Matrix<double, 3, 3> expected_covariance;
    expected_covariance << 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0;
    EXPECT_TRUE(expected_covariance.isApprox(math::rotate_point_covariance(covariance, rotation)));
}

TEST(rotate_point_covariance, yaw_180_x_covariance) {
    Eigen::Matrix<double, 3, 3> covariance;
    covariance << 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    Eigen::AngleAxisd rotation = Eigen::AngleAxisd(M_PI, Eigen::Vector3d::UnitZ());
    Eigen::Matrix<double, 3, 3> expected_covariance;
    expected_covariance << 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    EXPECT_TRUE(expected_covariance.isApprox(math::rotate_point_covariance(covariance, rotation)));
}

TEST(rotate_point_covariance, yaw_270_x_covariance) {
    Eigen::Matrix<double, 3, 3> covariance;
    covariance << 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    Eigen::AngleAxisd rotation = Eigen::AngleAxisd(1.5 * M_PI, Eigen::Vector3d::UnitZ());
    Eigen::Matrix<double, 3, 3> expected_covariance;
    expected_covariance << 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0;
    EXPECT_TRUE(expected_covariance.isApprox(math::rotate_point_covariance(covariance, rotation)));
}

TEST(rotate_point_covariance, yaw_45_x_covariance) {
    const double var = 0.5;
    Eigen::Matrix<double, 3, 3> covariance;
    covariance << var, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    Eigen::AngleAxisd rotation = Eigen::AngleAxisd(M_PI / 4.0, Eigen::Vector3d::UnitZ());
    Eigen::Vector3d xy = Eigen::Vector3d(1.0 / std::sqrt(2.0), 1.0 / std::sqrt(2.0), 0.0);
    Eigen::Matrix<double, 3, 3> expected_covariance = var * xy * xy.transpose();
    EXPECT_TRUE(expected_covariance.isApprox(math::rotate_point_covariance(covariance, rotation)));
}

TEST(transform_adjoint, rotation_only) {
    const Eigen::Isometry3d transform{Eigen::Quaterniond(0.17, 0.68, 0.55, 0.14).normalized()};
    Eigen::Matrix<double, 6, 6> transform_adjoint = math::transform_adjoint(transform);
    check_transform_adjoint_blocks(transform_adjoint, transform);
}

TEST(transform_adjoint, translation_only) {
    const Eigen::Isometry3d transform{Eigen::Translation<double, 3>{1.0, 2.0, 3.0}};
    const Eigen::Matrix<double, 6, 6> transform_adjoint = math::transform_adjoint(transform);
    check_transform_adjoint_blocks(transform_adjoint, transform);
}

TEST(transform_adjoint, transform) {
    const Eigen::Isometry3d transform =
            Eigen::Translation<double, 3>{1.0, 2.0, 3.0} * Eigen::Quaterniond(0.17, 0.68, 0.55, 0.14).normalized();
    const Eigen::Matrix<double, 6, 6> transform_adjoint = math::transform_adjoint(transform);
    check_transform_adjoint_blocks(transform_adjoint, transform);
}
