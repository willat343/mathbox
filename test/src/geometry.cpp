#include "mathbox/geometry.hpp"

#include <gtest/gtest.h>

#include <numbers>
#include <stdexcept>
#include <vector>

#include "mathbox/matrix_operations.hpp"

void check_adjoint_SE_rt_blocks(const Eigen::Matrix<double, 6, 6>& adjoint_SE, const Eigen::Isometry3d& transform) {
    const Eigen::Matrix3d R = transform.rotation();
    const Eigen::Vector3d t = transform.translation();
    const Eigen::Matrix3d t_SS = math::skew_symmetric_cross(t);
    const Eigen::Matrix3d t_SS_times_R = t_SS * R;
    const Eigen::Matrix3d top_left = adjoint_SE.block<3, 3>(0, 0);
    const Eigen::Matrix3d top_right = adjoint_SE.block<3, 3>(0, 3);
    const Eigen::Matrix3d bottom_left = adjoint_SE.block<3, 3>(3, 0);
    const Eigen::Matrix3d bottom_right = adjoint_SE.block<3, 3>(3, 3);
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
    EXPECT_TRUE(covariance.isApprox(math::change_tf_covariance_frame_tr(covariance, transform)));
    EXPECT_TRUE(covariance.isApprox(math::change_tf_covariance_frame_rt(covariance, transform)));
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
    const double dt = std::numbers::pi /
                      rotation_angle_rate;  // Too large a dt will invalidate test because the rotation will wrap
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

TEST(frame_from_directions, z_along_primary_and_x_along_secondary) {
    const Eigen::Vector3d primary{0.4, -0.2, 2.0};
    const Eigen::Vector3d secondary{3.0, 1.0, 0.5};
    const Eigen::Matrix3d M = math::frame_from_directions(primary, secondary);
    EXPECT_TRUE((M.transpose() * M).isIdentity(1.0e-12));
    EXPECT_NEAR(M.determinant(), 1.0, 1.0e-12);
    EXPECT_TRUE(M.col(2).isApprox(primary.normalized(), 1.0e-12));
    // The secondary direction is in the xz half-plane of positive x
    const Eigen::Vector3d secondary_M = M.transpose() * secondary;
    EXPECT_NEAR(secondary_M.y(), 0.0, 1.0e-12);
    EXPECT_GT(secondary_M.x(), 0.0);
}

TEST(frame_from_directions, identity_for_z_and_x) {
    EXPECT_TRUE(
            math::frame_from_directions(Eigen::Vector3d{2.0 * Eigen::Vector3d::UnitZ()}, Eigen::Vector3d{1.0, 0.0, 5.0})
                    .isIdentity(1.0e-12));
}

TEST(frame_from_directions, throws_for_zero_or_parallel_directions) {
    const Eigen::Vector3d x = Eigen::Vector3d::UnitX();
    const Eigen::Vector3d z = Eigen::Vector3d::UnitZ();
    const Eigen::Vector3d zero = Eigen::Vector3d::Zero();
    EXPECT_THROW(math::frame_from_directions(zero, x), std::invalid_argument);
    EXPECT_THROW(math::frame_from_directions(z, zero), std::invalid_argument);
    EXPECT_THROW(math::frame_from_directions(z, Eigen::Vector3d{-2.0 * z}), std::invalid_argument);
    EXPECT_THROW(math::frame_from_directions(z, z), std::invalid_argument);
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
    Eigen::AngleAxisd rotation = Eigen::AngleAxisd(std::numbers::pi / 2.0, Eigen::Vector3d::UnitZ());
    Eigen::Matrix<double, 3, 3> expected_covariance;
    expected_covariance << 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0;
    EXPECT_TRUE(expected_covariance.isApprox(math::rotate_point_covariance(covariance, rotation)));
}

TEST(rotate_point_covariance, yaw_180_x_covariance) {
    Eigen::Matrix<double, 3, 3> covariance;
    covariance << 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    Eigen::AngleAxisd rotation = Eigen::AngleAxisd(std::numbers::pi, Eigen::Vector3d::UnitZ());
    Eigen::Matrix<double, 3, 3> expected_covariance;
    expected_covariance << 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    EXPECT_TRUE(expected_covariance.isApprox(math::rotate_point_covariance(covariance, rotation)));
}

TEST(rotate_point_covariance, yaw_270_x_covariance) {
    Eigen::Matrix<double, 3, 3> covariance;
    covariance << 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    Eigen::AngleAxisd rotation = Eigen::AngleAxisd(1.5 * std::numbers::pi, Eigen::Vector3d::UnitZ());
    Eigen::Matrix<double, 3, 3> expected_covariance;
    expected_covariance << 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0;
    EXPECT_TRUE(expected_covariance.isApprox(math::rotate_point_covariance(covariance, rotation)));
}

TEST(rotate_point_covariance, yaw_45_x_covariance) {
    const double var = 0.5;
    Eigen::Matrix<double, 3, 3> covariance;
    covariance << var, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    Eigen::AngleAxisd rotation = Eigen::AngleAxisd(std::numbers::pi / 4.0, Eigen::Vector3d::UnitZ());
    Eigen::Vector3d xy = Eigen::Vector3d(1.0 / std::sqrt(2.0), 1.0 / std::sqrt(2.0), 0.0);
    Eigen::Matrix<double, 3, 3> expected_covariance = var * xy * xy.transpose();
    EXPECT_TRUE(expected_covariance.isApprox(math::rotate_point_covariance(covariance, rotation)));
}

TEST(adjoint_SE, rotation_only) {
    const Eigen::Isometry3d transform{Eigen::Quaterniond(0.17, 0.68, 0.55, 0.14).normalized()};
    Eigen::Matrix<double, 6, 6> adjoint_SE = math::adjoint_SE_rt(transform);
    check_adjoint_SE_rt_blocks(adjoint_SE, transform);
}

TEST(adjoint_SE, translation_only) {
    const Eigen::Isometry3d transform{Eigen::Translation<double, 3>{1.0, 2.0, 3.0}};
    const Eigen::Matrix<double, 6, 6> adjoint_SE = math::adjoint_SE_rt(transform);
    check_adjoint_SE_rt_blocks(adjoint_SE, transform);
}

TEST(adjoint_SE, transform) {
    const Eigen::Isometry3d transform =
            Eigen::Translation<double, 3>{1.0, 2.0, 3.0} * Eigen::Quaterniond(0.17, 0.68, 0.55, 0.14).normalized();
    const Eigen::Matrix<double, 6, 6> adjoint_SE = math::adjoint_SE_rt(transform);
    check_adjoint_SE_rt_blocks(adjoint_SE, transform);
}

TEST(angle_between, parallel_3d) {
    const Eigen::Vector3d u(1.0, 2.0, 3.0);
    const Eigen::Vector3d v = 2.0 * u;
    EXPECT_NEAR(math::angle_between(u, v), 0.0, 1.0e-12);
}

TEST(angle_between, opposite_3d) {
    const Eigen::Vector3d u(1.0, 2.0, 3.0);
    const Eigen::Vector3d v = -u;
    EXPECT_NEAR(math::angle_between(u, v), std::numbers::pi, 1.0e-12);
}

TEST(angle_between, orthogonal_3d) {
    const Eigen::Vector3d u = Eigen::Vector3d::UnitX();
    const Eigen::Vector3d v = Eigen::Vector3d::UnitY();
    EXPECT_NEAR(math::angle_between(u, v), std::numbers::pi / 2.0, 1.0e-12);
}

TEST(angle_between, orthogonal_2d) {
    const Eigen::Vector2d u(1.0, 0.0);
    const Eigen::Vector2d v(0.0, 1.0);
    EXPECT_NEAR(math::angle_between(u, v), std::numbers::pi / 2.0, 1.0e-12);
}

TEST(deg2rad, 90_degrees) {
    EXPECT_NEAR(math::deg2rad(90.0), std::numbers::pi / 2.0, 1.0e-12);
}

TEST(deg2rad, 180_degrees) {
    EXPECT_NEAR(math::deg2rad(180.0), std::numbers::pi, 1.0e-12);
}

TEST(rad2deg, pi_over_2) {
    EXPECT_NEAR(math::rad2deg(std::numbers::pi / 2.0), 90.0, 1.0e-12);
}

TEST(rad2deg, pi) {
    EXPECT_NEAR(math::rad2deg(std::numbers::pi), 180.0, 1.0e-12);
}

TEST(deg2rad, round_trip_with_rad2deg) {
    const double degrees = 37.5;
    EXPECT_NEAR(math::rad2deg(math::deg2rad(degrees)), degrees, 1.0e-12);
}

TEST(rpy, identity) {
    const Eigen::Quaterniond q = Eigen::Quaterniond::Identity();
    EXPECT_TRUE(math::rpy(q).isApprox(Eigen::Vector3d::Zero(), 1.0e-12));
}

TEST(rpy, roll_only) {
    const double roll = 0.4;
    const Eigen::Quaterniond q(Eigen::AngleAxisd(roll, Eigen::Vector3d::UnitX()));
    const Eigen::Vector3d angles = math::rpy(q);
    EXPECT_NEAR(angles[0], roll, 1.0e-9);
    EXPECT_NEAR(angles[1], 0.0, 1.0e-9);
    EXPECT_NEAR(angles[2], 0.0, 1.0e-9);
}

TEST(rpy, yaw_only) {
    const double yaw = 0.3;
    const Eigen::Quaterniond q(Eigen::AngleAxisd(yaw, Eigen::Vector3d::UnitZ()));
    const Eigen::Vector3d angles = math::rpy(q);
    EXPECT_NEAR(angles[0], 0.0, 1.0e-9);
    EXPECT_NEAR(angles[1], 0.0, 1.0e-9);
    EXPECT_NEAR(angles[2], yaw, 1.0e-9);
}

TEST(rotation_from_directions, recovers_rotation_from_consistent_pairs) {
    const std::vector<Eigen::Quaterniond> rotations{Eigen::Quaterniond::Identity(),
            Eigen::Quaterniond{Eigen::AngleAxisd{std::numbers::pi, Eigen::Vector3d::UnitX()}},
            Eigen::Quaterniond{0.8, 0.1, 0.05, 0.2}.normalized(),
            Eigen::Quaterniond{-0.1, 0.7, -0.6, 0.3}.normalized()};
    const Eigen::Vector3d primary_B{0.1, -0.2, 3.0};
    const Eigen::Vector3d secondary_B{2.0, 0.5, -0.4};
    for (const Eigen::Quaterniond& rotation : rotations) {
        SCOPED_TRACE(rotation.coeffs().transpose());
        const Eigen::Matrix3d R_A_B = rotation.toRotationMatrix();
        // The lengths of the directions do not matter
        const Eigen::Vector3d primary_A = 0.5 * R_A_B * primary_B;
        const Eigen::Vector3d secondary_A = 4.0 * R_A_B * secondary_B;
        EXPECT_TRUE(math::rotation_from_directions(primary_A, secondary_A, primary_B, secondary_B)
                        .isApprox(R_A_B, 1.0e-12));
    }
}

TEST(rotation_from_directions, primary_exact_and_secondary_about_it_for_inconsistent_pairs) {
    // The angle between the directions differs between the frames
    const Eigen::Vector3d primary_A = Eigen::Vector3d::UnitZ();
    const Eigen::Vector3d secondary_A{1.0, 0.0, 0.3};
    const Eigen::Vector3d primary_B = Eigen::Vector3d{0.2, -0.1, 1.0}.normalized();
    const Eigen::Vector3d secondary_B{0.0, 1.0, -0.5};
    const Eigen::Matrix3d R_A_B = math::rotation_from_directions(primary_A, secondary_A, primary_B, secondary_B);
    EXPECT_TRUE((R_A_B.transpose() * R_A_B).isIdentity(1.0e-12));
    EXPECT_NEAR(R_A_B.determinant(), 1.0, 1.0e-12);
    EXPECT_TRUE((R_A_B * primary_B).isApprox(primary_A, 1.0e-12));
    // The secondary direction is mapped into the half-plane of the primary and secondary directions of A
    const Eigen::Vector3d secondary = R_A_B * secondary_B;
    EXPECT_NEAR(secondary.y(), 0.0, 1.0e-12);
    EXPECT_GT(secondary.x(), 0.0);
}

TEST(rotation_from_directions, opposite_primary_directions_keep_secondary_heading) {
    // A primary direction close to z in B and -z in A (e.g. gravity measured in a body frame whose z axis is down,
    // in a z-up frame): the x axis of B keeps its heading in A, whatever the direction of the small tilt
    const std::vector<Eigen::Vector3d> primaries_B{Eigen::Vector3d::UnitZ(), Eigen::Vector3d{0.01, 0.02, 1.0},
            Eigen::Vector3d{-0.02, 0.01, 1.0}, Eigen::Vector3d{0.03, -0.01, 1.0}};
    for (const Eigen::Vector3d& primary_B : primaries_B) {
        SCOPED_TRACE(primary_B.transpose());
        const Eigen::Matrix3d R_A_B = math::rotation_from_directions(Eigen::Vector3d{-Eigen::Vector3d::UnitZ()},
                Eigen::Vector3d{Eigen::Vector3d::UnitX()}, primary_B, Eigen::Vector3d{Eigen::Vector3d::UnitX()});
        EXPECT_TRUE((R_A_B * primary_B.normalized()).isApprox(-Eigen::Vector3d::UnitZ(), 1.0e-12));
        EXPECT_NEAR(R_A_B(1, 0), 0.0, 1.0e-12);
        EXPECT_GT(R_A_B(0, 0), 0.99);
        EXPECT_LT(R_A_B(2, 2), -0.99);
    }
}

TEST(rotation_from_directions, throws_for_zero_or_parallel_directions) {
    const Eigen::Vector3d x = Eigen::Vector3d::UnitX();
    const Eigen::Vector3d z = Eigen::Vector3d::UnitZ();
    EXPECT_THROW(math::rotation_from_directions(Eigen::Vector3d{Eigen::Vector3d::Zero()}, x, z, x),
            std::invalid_argument);
    EXPECT_THROW(math::rotation_from_directions(z, x, z, Eigen::Vector3d{-2.0 * z}), std::invalid_argument);
}

TEST(quaternion_to_rotation_vector, matches_angle_axis) {
    for (const double angle : {0.3, 1.0, 2.5, 3.1}) {
        SCOPED_TRACE(angle);
        const Eigen::Vector3d axis = Eigen::Vector3d(0.4, -1.0, 0.7).normalized();
        EXPECT_TRUE(math::quaternion_to_rotation_vector(Eigen::Quaterniond(Eigen::AngleAxisd(angle, axis)))
                        .isApprox(angle * axis));
    }
}

TEST(quaternion_to_rotation_vector, shortest_rotation_for_either_sign) {
    const Eigen::Vector3d axis = Eigen::Vector3d(0.4, -1.0, 0.7).normalized();
    const Eigen::Quaterniond rotation(Eigen::AngleAxisd(0.8, axis));
    EXPECT_TRUE(math::quaternion_to_rotation_vector(Eigen::Quaterniond(-rotation.coeffs())).isApprox(0.8 * axis));
    // A rotation by 4 rad is the shortest rotation by 2 pi - 4 rad about the opposite axis
    EXPECT_TRUE(math::quaternion_to_rotation_vector(Eigen::Quaterniond(Eigen::AngleAxisd(4.0, axis)))
                    .isApprox((4.0 - 2.0 * std::numbers::pi) * axis));
}

TEST(quaternion_to_rotation_vector, small_angle) {
    const Eigen::Vector3d axis = Eigen::Vector3d(0.4, -1.0, 0.7).normalized();
    const Eigen::Quaterniond rotation(Eigen::AngleAxisd(1.0e-9, axis));
    EXPECT_TRUE(math::quaternion_to_rotation_vector(rotation).isApprox(1.0e-9 * axis));
    EXPECT_TRUE(math::quaternion_to_rotation_vector(Eigen::Quaterniond(-rotation.coeffs())).isApprox(1.0e-9 * axis));
    EXPECT_TRUE(math::quaternion_to_rotation_vector(Eigen::Quaterniond::Identity()).isZero(0.0));
}

TEST(rotation_vector_to_quaternion, matches_angle_axis) {
    const Eigen::Vector3d axis = Eigen::Vector3d(0.4, -1.0, 0.7).normalized();
    for (const double angle : {0.3, 1.0, 2.5, 3.1}) {
        SCOPED_TRACE(angle);
        EXPECT_TRUE(math::rotation_vector_to_quaternion<double>(angle * axis)
                        .coeffs()
                        .isApprox(Eigen::Quaterniond(Eigen::AngleAxisd(angle, axis)).coeffs()));
    }
}

TEST(rotation_vector_to_quaternion, small_angle) {
    const Eigen::Vector3d axis = Eigen::Vector3d(0.4, -1.0, 0.7).normalized();
    EXPECT_TRUE(math::rotation_vector_to_quaternion<double>(1.0e-9 * axis)
                    .coeffs()
                    .isApprox(Eigen::Quaterniond(Eigen::AngleAxisd(1.0e-9, axis)).coeffs()));
    EXPECT_TRUE(math::rotation_vector_to_quaternion<double>(Eigen::Vector3d::Zero())
                    .coeffs()
                    .isApprox(Eigen::Quaterniond::Identity().coeffs()));
}

TEST(rotation_vector_to_quaternion, inverts_quaternion_to_rotation_vector) {
    const Eigen::Vector3d axis = Eigen::Vector3d(0.4, -1.0, 0.7).normalized();
    for (const double angle : {1.0e-9, 0.3, 3.1}) {
        SCOPED_TRACE(angle);
        EXPECT_TRUE(math::quaternion_to_rotation_vector(math::rotation_vector_to_quaternion<double>(angle * axis))
                        .isApprox(angle * axis));
    }
}

TEST(so_cross, so2) {
    const Eigen::Matrix<double, 1, 1> w(2.0);
    const Eigen::Vector2d v(3.0, 4.0);
    const Eigen::Vector2d expected(-2.0 * 4.0, 2.0 * 3.0);
    EXPECT_TRUE(math::so_cross(w, v).isApprox(expected));
}

TEST(so_cross, so3_matches_eigen_cross) {
    const Eigen::Vector3d w(1.0, 2.0, 3.0);
    const Eigen::Vector3d v(4.0, 5.0, 6.0);
    EXPECT_TRUE(math::so_cross(w, v).isApprox(w.cross(v)));
}

// Regression test: so_skew's 2D specialization used to build a matrix that was not even skew-symmetric, so
// so_skew(w) * v did not equal so_cross(w, v).
TEST(so_skew, so2_matches_so_cross) {
    const Eigen::Matrix<double, 1, 1> w(2.0);
    const Eigen::Vector2d v(3.0, 4.0);
    const Eigen::Matrix2d skew_w = math::so_skew(w);
    EXPECT_TRUE((skew_w * v).isApprox(math::so_cross(w, v)));
}

TEST(so_skew, so2_matrix_values) {
    const Eigen::Matrix<double, 1, 1> w(2.0);
    Eigen::Matrix2d expected;
    expected << 0.0, -2.0, 2.0, 0.0;
    EXPECT_TRUE(math::so_skew(w).isApprox(expected));
}

TEST(so_skew, so3_matches_skew_symmetric_cross) {
    const Eigen::Vector3d w(1.0, 2.0, 3.0);
    EXPECT_TRUE(math::so_skew(w).isApprox(math::skew_symmetric_cross(w)));
}

// Regression test: so_from_skew's 2D specialization used to read the wrong matrix entry (self-consistent with the
// buggy so_skew above, but not the true so(2) vee map).
TEST(so_from_skew, so2_round_trip) {
    const Eigen::Matrix<double, 1, 1> w(2.0);
    const Eigen::Matrix2d skew_w = math::so_skew(w);
    EXPECT_TRUE(math::so_from_skew(skew_w).isApprox(w));
}

TEST(so_from_skew, so3_round_trip) {
    const Eigen::Vector3d w(1.0, 2.0, 3.0);
    const Eigen::Matrix3d skew_w = math::so_skew(w);
    EXPECT_TRUE(math::so_from_skew(skew_w).isApprox(w));
}

// Regression test: to_pose_2D used to throw when the axis matched (the valid case) and stay silent when it didn't
// (the actual error case), due to a missing negation.
TEST(to_pose_2D, matching_axis_does_not_throw) {
    const Eigen::Vector3d translation(1.0, 2.0, 3.0);
    const double angle = 0.7;
    const math::Pose<3> pose = Eigen::Translation3d(translation) * Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitZ());
    const math::Pose<2> pose_2d = math::to_pose_2D(pose, Eigen::Vector3d::UnitZ());
    EXPECT_NEAR(pose_2d.translation().x(), translation.x(), 1.0e-12);
    EXPECT_NEAR(pose_2d.translation().y(), translation.y(), 1.0e-12);
    EXPECT_NEAR(Eigen::Rotation2Dd(pose_2d.rotation()).angle(), angle, 1.0e-12);
}

TEST(to_pose_2D, mismatched_axis_throws) {
    const math::Pose<3> pose{Eigen::AngleAxisd(0.7, Eigen::Vector3d::UnitZ())};
    EXPECT_THROW(math::to_pose_2D(pose, Eigen::Vector3d::UnitX()), std::invalid_argument);
}

TEST(to_pose_2D, default_axis_skips_check) {
    const math::Pose<3> pose{Eigen::AngleAxisd(0.7, Eigen::Vector3d::UnitX())};
    EXPECT_NO_THROW(math::to_pose_2D(pose));
}

TEST(to_pose_ND, D3_is_passthrough) {
    const math::Pose<3> pose = Eigen::Translation3d(1.0, 2.0, 3.0) * Eigen::AngleAxisd(0.5, Eigen::Vector3d::UnitY());
    EXPECT_TRUE(math::to_pose_ND<3>(pose).isApprox(pose));
}

TEST(to_pose_ND, D2_matches_to_pose_2D) {
    const Eigen::Vector3d translation(1.0, 2.0, 3.0);
    const double angle = 0.7;
    const math::Pose<3> pose = Eigen::Translation3d(translation) * Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitZ());
    const math::Pose<2> pose_nd = math::to_pose_ND<2>(pose, Eigen::Vector3d::UnitZ());
    const math::Pose<2> pose_2d = math::to_pose_2D(pose, Eigen::Vector3d::UnitZ());
    EXPECT_TRUE(pose_nd.isApprox(pose_2d));
}

/**
 * @brief Apply a perturbation to a pose, to first order in the position tangent vector.
 */
Eigen::Isometry3d perturb_pose(const Eigen::Isometry3d& pose, const math::PosePerturbation& perturbation,
        const Eigen::Vector<double, 6>& tangent) {
    const Eigen::Vector3d rho = tangent.head<3>();
    const Eigen::Vector3d theta = tangent.tail<3>();
    const Eigen::Matrix3d exp_theta = Eigen::AngleAxisd(theta.norm(), theta.normalized()).toRotationMatrix();
    const Eigen::Matrix3d R = pose.rotation();
    const Eigen::Vector3d p = pose.translation();
    Eigen::Isometry3d perturbed = Eigen::Isometry3d::Identity();
    perturbed.linear() = perturbation.side == math::ApplicationSide::LEFT ? exp_theta * R : R * exp_theta;
    if (perturbation.type == math::SEType::SPLIT) {
        perturbed.translation() = p + rho;
    } else if (perturbation.side == math::ApplicationSide::LEFT) {
        perturbed.translation() = exp_theta * p + rho;
    } else {
        perturbed.translation() = p + R * rho;
    }
    return perturbed;
}

std::vector<math::PosePerturbation> pose_perturbations() {
    return {{math::SEType::JOINT, math::ApplicationSide::LEFT}, {math::SEType::JOINT, math::ApplicationSide::RIGHT},
            {math::SEType::SPLIT, math::ApplicationSide::LEFT}, {math::SEType::SPLIT, math::ApplicationSide::RIGHT}};
}

Eigen::Isometry3d perturbation_test_pose() {
    Eigen::Isometry3d pose = Eigen::Isometry3d::Identity();
    pose.linear() = Eigen::AngleAxisd(0.7, Eigen::Vector3d(1.0, -2.0, 3.0).normalized()).toRotationMatrix();
    pose.translation() = Eigen::Vector3d(4.0, -5.0, 6.0);
    return pose;
}

TEST(PosePerturbation, equality) {
    const math::PosePerturbation perturbation{math::SEType::SPLIT, math::ApplicationSide::RIGHT};
    EXPECT_TRUE(perturbation == (math::PosePerturbation{math::SEType::SPLIT, math::ApplicationSide::RIGHT}));
    EXPECT_FALSE(perturbation == (math::PosePerturbation{math::SEType::SPLIT, math::ApplicationSide::LEFT}));
    EXPECT_FALSE(perturbation == (math::PosePerturbation{math::SEType::JOINT, math::ApplicationSide::RIGHT}));
}

TEST(pose_perturbation_jacobian, same_perturbation_is_identity) {
    for (const math::PosePerturbation& perturbation : pose_perturbations()) {
        EXPECT_TRUE(math::pose_perturbation_jacobian(perturbation_test_pose(), perturbation, perturbation)
                        .isIdentity(1.0e-12));
    }
}

TEST(pose_perturbation_jacobian, inverse) {
    for (const math::PosePerturbation& perturbation : pose_perturbations()) {
        for (const math::PosePerturbation& with_respect_to : pose_perturbations()) {
            EXPECT_TRUE((math::pose_perturbation_jacobian(perturbation_test_pose(), perturbation, with_respect_to) *
                         math::pose_perturbation_jacobian(perturbation_test_pose(), with_respect_to, perturbation))
                            .isIdentity(1.0e-12));
        }
    }
}

TEST(pose_perturbation_jacobian, joint_left_with_respect_to_joint_right_is_adjoint) {
    EXPECT_TRUE(math::pose_perturbation_jacobian(perturbation_test_pose(),
            {math::SEType::JOINT, math::ApplicationSide::LEFT}, {math::SEType::JOINT, math::ApplicationSide::RIGHT})
                    .isApprox(math::adjoint_SE_tr(perturbation_test_pose())));
}

TEST(pose_perturbation_jacobian, perturbed_poses_match) {
    const Eigen::Isometry3d pose = perturbation_test_pose();
    const Eigen::Vector<double, 6> tangent =
            1.0e-6 * (Eigen::Vector<double, 6>() << 1.0, -2.0, 3.0, -1.5, 2.5, 0.5).finished();
    for (const math::PosePerturbation& perturbation : pose_perturbations()) {
        for (const math::PosePerturbation& with_respect_to : pose_perturbations()) {
            const Eigen::Vector<double, 6> converted_tangent =
                    math::pose_perturbation_jacobian(pose, perturbation, with_respect_to) * tangent;
            // Equal to first order, so the difference is of the order of the squared tangent vector
            EXPECT_TRUE(perturb_pose(pose, perturbation, converted_tangent)
                            .matrix()
                            .isApprox(perturb_pose(pose, with_respect_to, tangent).matrix(), 1.0e-10));
        }
    }
}

TEST(pose_perturbation_jacobian, planar) {
    Eigen::Isometry2d pose = Eigen::Isometry2d::Identity();
    pose.linear() = Eigen::Rotation2Dd(0.7).toRotationMatrix();
    pose.translation() = Eigen::Vector2d(4.0, -5.0);
    // The orientation tangent of a planar pose is the same on both sides
    EXPECT_TRUE(math::pose_perturbation_jacobian(pose, {math::SEType::SPLIT, math::ApplicationSide::LEFT},
            {math::SEType::SPLIT, math::ApplicationSide::RIGHT})
                    .isIdentity(1.0e-12));
    EXPECT_TRUE(math::pose_perturbation_jacobian(pose, {math::SEType::JOINT, math::ApplicationSide::LEFT},
            {math::SEType::JOINT, math::ApplicationSide::RIGHT})
                    .isApprox(math::adjoint_SE_tr(pose)));
}
