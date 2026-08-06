#include "mathbox/matrix_properties.hpp"

#include <gtest/gtest.h>

#include <cmath>

TEST(has_positive_diagonals, all_positive) {
    Eigen::Matrix2d m;
    m << 1.0, -5.0, -5.0, 2.0;
    EXPECT_TRUE(math::has_positive_diagonals(m));
}

TEST(has_positive_diagonals, zero_diagonal_is_not_positive) {
    Eigen::Matrix2d m;
    m << 1.0, 0.0, 0.0, 0.0;
    EXPECT_FALSE(math::has_positive_diagonals(m));
}

TEST(is_symmetric, symmetric_matrix) {
    Eigen::Matrix2d m;
    m << 1.0, 2.0, 2.0, 3.0;
    EXPECT_TRUE(math::is_symmetric(m));
}

TEST(is_symmetric, asymmetric_matrix) {
    Eigen::Matrix2d m;
    m << 1.0, 2.0, 4.0, 3.0;
    EXPECT_FALSE(math::is_symmetric(m));
}

TEST(is_symmetric, asymmetric_within_precision) {
    Eigen::Matrix2d m;
    m << 1.0, 2.0, 2.0 + 1.0e-9, 3.0;
    EXPECT_FALSE(math::is_symmetric(m));
    EXPECT_TRUE(math::is_symmetric(m, 1.0e-6));
}

TEST(is_skew_symmetric, skew_symmetric_matrix) {
    Eigen::Matrix2d m;
    m << 0.0, -2.0, 2.0, 0.0;
    EXPECT_TRUE(math::is_skew_symmetric(m));
}

TEST(is_skew_symmetric, symmetric_matrix_is_not_skew_symmetric) {
    Eigen::Matrix2d m;
    m << 1.0, 2.0, 2.0, 3.0;
    EXPECT_FALSE(math::is_skew_symmetric(m));
}

TEST(is_positive_definite, identity) {
    EXPECT_TRUE(math::is_positive_definite(Eigen::Matrix2d(Eigen::Matrix2d::Identity())));
}

TEST(is_positive_definite, zero_matrix_is_not_positive_definite) {
    EXPECT_FALSE(math::is_positive_definite(Eigen::Matrix2d(Eigen::Matrix2d::Zero())));
}

TEST(is_positive_definite, asymmetric_matrix_is_not_positive_definite) {
    Eigen::Matrix2d m;
    m << 1.0, 2.0, 4.0, 3.0;
    EXPECT_FALSE(math::is_positive_definite(m));
}

TEST(is_positive_semidefinite, identity) {
    EXPECT_TRUE(math::is_positive_semidefinite(Eigen::Matrix2d(Eigen::Matrix2d::Identity())));
}

TEST(is_positive_semidefinite, zero_matrix_is_positive_semidefinite) {
    EXPECT_TRUE(math::is_positive_semidefinite(Eigen::Matrix2d(Eigen::Matrix2d::Zero())));
}

TEST(is_positive_semidefinite, negative_definite_matrix) {
    EXPECT_FALSE(math::is_positive_semidefinite(Eigen::Matrix2d(-Eigen::Matrix2d::Identity())));
}

TEST(is_upper_triangular, upper_triangular_matrix) {
    Eigen::Matrix2d m;
    m << 1.0, 2.0, 0.0, 3.0;
    EXPECT_TRUE(math::is_upper_triangular(m));
}

// Regression test: the loop bound used to skip the last row entirely, so a non-zero entry below the diagonal in the
// final row was never checked.
TEST(is_upper_triangular, nonzero_in_last_row_is_not_upper_triangular) {
    Eigen::Matrix2d m;
    m << 1.0, 0.0, 5.0, 1.0;
    EXPECT_FALSE(math::is_upper_triangular(m));
}

TEST(is_upper_triangular, nonzero_in_last_row_3x3_is_not_upper_triangular) {
    Eigen::Matrix3d m;
    m << 1.0, 2.0, 3.0, 0.0, 4.0, 5.0, 0.0, 7.0, 6.0;
    EXPECT_FALSE(math::is_upper_triangular(m));
}

TEST(is_upper_triangular, within_precision) {
    Eigen::Matrix2d m;
    m << 1.0, 2.0, 1.0e-9, 3.0;
    EXPECT_FALSE(math::is_upper_triangular(m));
    EXPECT_TRUE(math::is_upper_triangular(m, 1.0e-6));
}

TEST(num_zero_columns, counts_zero_columns) {
    Eigen::Matrix3d m;
    m << 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 3.0, 0.0, 0.0;
    EXPECT_EQ(math::num_zero_columns(m), 2);
}

TEST(num_zero_columns, no_zero_columns) {
    EXPECT_EQ(math::num_zero_columns(Eigen::Matrix3d(Eigen::Matrix3d::Identity())), 0);
}

TEST(num_zero_rows, counts_zero_rows) {
    Eigen::Matrix3d m;
    m << 1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    EXPECT_EQ(math::num_zero_rows(m), 2);
}

TEST(num_zero_rows, no_zero_rows) {
    EXPECT_EQ(math::num_zero_rows(Eigen::Matrix3d(Eigen::Matrix3d::Identity())), 0);
}

TEST(relative_asymmetry, symmetric_matrix_is_zero) {
    Eigen::Matrix2d m;
    m << 1.0, 2.0, 2.0, 3.0;
    EXPECT_DOUBLE_EQ(math::relative_asymmetry(m), 0.0);
}

// Regression test: the open-coded form of this check divided by the matrix norm unguarded, so an exactly-zero matrix
// produced NaN (and any comparison against NaN is false, silently bypassing the symmetry check).
TEST(relative_asymmetry, zero_matrix_is_zero_not_nan) {
    const double asymmetry = math::relative_asymmetry(Eigen::Matrix2d(Eigen::Matrix2d::Zero()));
    EXPECT_FALSE(std::isnan(asymmetry));
    EXPECT_DOUBLE_EQ(asymmetry, 0.0);
}

TEST(relative_asymmetry, asymmetric_matrix_is_positive) {
    Eigen::Matrix2d m;
    m << 1.0, 2.0, 4.0, 3.0;
    // ||m - m^T||_F = ||[[0, -2], [2, 0]]||_F = sqrt(8), ||m||_F = sqrt(30).
    EXPECT_NEAR(math::relative_asymmetry(m), std::sqrt(8.0) / std::sqrt(30.0), 1.0e-12);
}

TEST(relative_asymmetry, scale_invariant) {
    Eigen::Matrix2d m;
    m << 1.0, 2.0, 4.0, 3.0;
    EXPECT_NEAR(math::relative_asymmetry(m), math::relative_asymmetry(Eigen::Matrix2d(1000.0 * m)), 1.0e-12);
}
