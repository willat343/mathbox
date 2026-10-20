#include "mathbox/matrix_operations.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>

TEST(cumulative_row_left_sum, 2x2) {
    Eigen::Matrix2d m;
    m << 1.0, 2.0, 3.0, 4.0;
    Eigen::Matrix2d expected;
    expected << 3.0, 2.0, 7.0, 4.0;
    const Eigen::Matrix2d result = math::cumulative_row_left_sum(m);
    EXPECT_TRUE(result.isApprox(expected));
}

TEST(cumulative_col_top_sum, 2x2) {
    Eigen::Matrix2d m;
    m << 1.0, 2.0, 3.0, 4.0;
    Eigen::Matrix2d expected;
    expected << 4.0, 6.0, 3.0, 4.0;
    const Eigen::Matrix2d result = math::cumulative_col_top_sum(m);
    EXPECT_TRUE(result.isApprox(expected));
}

TEST(make_symmetric, 2x2) {
    Eigen::Matrix2d m;
    m << 1.0, 2.0, 4.0, 3.0;
    Eigen::Matrix2d expected;
    expected << 1.0, 3.0, 3.0, 3.0;
    const Eigen::Matrix2d result = math::make_symmetric(m);
    EXPECT_TRUE(result.isApprox(expected));
    EXPECT_TRUE(result.isApprox(result.transpose()));
}

TEST(remove_rows_by_threshold, keeps_rows_greater_equal_threshold) {
    Eigen::MatrixXd m(4, 2);
    m << 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0;
    Eigen::VectorXd v(4);
    v << 0.0, 1.0, 2.0, 3.0;
    const Eigen::MatrixXd result = math::remove_rows_by_threshold(m, v, 2.0);
    Eigen::MatrixXd expected(2, 2);
    expected << 3.0, 3.0, 4.0, 4.0;
    EXPECT_TRUE(result.isApprox(expected));
}

TEST(skew_symmetric_cross, cross_product_equivalence) {
    const Eigen::Vector3d x(1.0, 2.0, 3.0);
    const Eigen::Vector3d y(4.0, 5.0, 6.0);
    const Eigen::Matrix3d skew_x = math::skew_symmetric_cross(x);
    EXPECT_TRUE((skew_x * y).isApprox(x.cross(y)));
}

TEST(skew_symmetric_cross, matrix_values) {
    const Eigen::Vector3d x(1.0, 2.0, 3.0);
    Eigen::Matrix3d expected;
    expected << 0.0, -3.0, 2.0, 3.0, 0.0, -1.0, -2.0, 1.0, 0.0;
    const Eigen::Matrix3d skew_x = math::skew_symmetric_cross(x);
    EXPECT_TRUE(skew_x.isApprox(expected));
}

TEST(schur_complement, three_variables_marginalise_one) {
    Eigen::MatrixXd H(3, 3);
    H << 4.0, 1.0, 2.0, 1.0, 5.0, 1.0, 2.0, 1.0, 6.0;
    Eigen::VectorXd b(3);
    b << 1.0, 2.0, 3.0;
    Eigen::MatrixXd H_p;
    Eigen::VectorXd b_p;
    const double damping = math::schur_complement(H, b, 1, H_p, b_p, 0.0, 1.0e-9, false);

    Eigen::MatrixXd H_p_expected(2, 2);
    H_p_expected << 4.75, 0.5, 0.5, 5.0;
    Eigen::VectorXd b_p_expected(2);
    b_p_expected << 1.75, 2.5;

    EXPECT_TRUE(H_p.isApprox(H_p_expected));
    EXPECT_TRUE(b_p.isApprox(b_p_expected));
    EXPECT_DOUBLE_EQ(damping, 0.0);
}

TEST(schur_complement, marginalise_leaving_single_variable) {
    // Regression test: upper_block_size == b.size() - 1 (lower_block_size == 1) used to violate a debug-only assert.
    Eigen::MatrixXd H(3, 3);
    H << 4.0, 1.0, 2.0, 1.0, 5.0, 1.0, 2.0, 1.0, 6.0;
    Eigen::VectorXd b(3);
    b << 1.0, 2.0, 3.0;
    Eigen::MatrixXd H_p;
    Eigen::VectorXd b_p;
    const double damping = math::schur_complement(H, b, 2, H_p, b_p, 0.0, 1.0e-9, false);

    EXPECT_NEAR(H_p(0, 0), 94.0 / 19.0, 1.0e-9);
    EXPECT_NEAR(b_p(0), 44.0 / 19.0, 1.0e-9);
    EXPECT_DOUBLE_EQ(damping, 0.0);
}

TEST(schur_complement, zero_h_p_does_not_produce_nan) {
    // Regression test: H_p.norm() == 0 used to cause a 0/0 division to NaN in the symmetry check.
    Eigen::MatrixXd H(2, 2);
    H << 4.0, 2.0, 2.0, 1.0;
    Eigen::VectorXd b(2);
    b << 1.0, 1.0;
    Eigen::MatrixXd H_p;
    Eigen::VectorXd b_p;
    EXPECT_NO_THROW(math::schur_complement(H, b, 1, H_p, b_p, 0.0, 1.0e-9, false));
    EXPECT_TRUE(H_p.isApprox(Eigen::MatrixXd::Zero(1, 1)));
    EXPECT_FALSE(std::isnan(H_p(0, 0)));
}

// Deterministic non-dyadic symmetric positive definite matrix, so that the scaled and unscaled paths round
// differently and can only be compared approximately.
Eigen::MatrixXd well_conditioned_spd(const int size) {
    Eigen::MatrixXd M(size, size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            M(i, j) = std::sin(static_cast<double>(size * i + j + 1));
        }
    }
    return M * M.transpose() + static_cast<double>(size) * Eigen::MatrixXd::Identity(size, size);
}

Eigen::VectorXd deterministic_vector(const int size) {
    Eigen::VectorXd v(size);
    for (int i = 0; i < size; ++i) {
        v[i] = std::cos(static_cast<double>(i + 1));
    }
    return v;
}

TEST(schur_complement, jacobi_scaling_is_noop_on_well_scaled_matrix) {
    // The diagonal of H_mm is 4, so the scaling is exactly 2 and both paths are exact in binary floating point.
    Eigen::MatrixXd H(3, 3);
    H << 4.0, 1.0, 2.0, 1.0, 5.0, 1.0, 2.0, 1.0, 6.0;
    Eigen::VectorXd b(3);
    b << 1.0, 2.0, 3.0;
    Eigen::MatrixXd H_p;
    Eigen::VectorXd b_p;
    const double damping = math::schur_complement(H, b, 1, H_p, b_p, 0.0, 1.0e-9, true);

    EXPECT_DOUBLE_EQ(H_p(0, 0), 4.75);
    EXPECT_DOUBLE_EQ(H_p(0, 1), 0.5);
    EXPECT_DOUBLE_EQ(H_p(1, 0), 0.5);
    EXPECT_DOUBLE_EQ(H_p(1, 1), 5.0);
    EXPECT_DOUBLE_EQ(b_p(0), 1.75);
    EXPECT_DOUBLE_EQ(b_p(1), 2.5);
    EXPECT_DOUBLE_EQ(damping, 0.0);
}

TEST(schur_complement, jacobi_scaling_badly_scaled_matrix_analytic) {
    // The diagonal of H_mm spans 2^30, and all entries are dyadic so that both paths are exact in binary floating
    // point. The two scalings differ (2^20 and 2^5), so omitting the D^-1 on b_m changes b_p by a factor of 2^19.
    const double p5 = 32.0, p6 = 64.0, p10 = 1024.0, p20 = 1048576.0, p21 = 2097152.0, p40 = p20 * p20;
    Eigen::MatrixXd H(4, 4);
    H << p40, 0.0, p20, p21, 0.0, p10, p5, p6, p20, p5, 8.0, 0.0, p21, p6, 0.0, 16.0;
    Eigen::VectorXd b(4);
    b << p21, p6, 5.0, 10.0;

    for (const bool jacobi_scaling : {false, true}) {
        Eigen::MatrixXd H_p;
        Eigen::VectorXd b_p;
        math::schur_complement(H, b, 2, H_p, b_p, 0.0, 1.0e-9, jacobi_scaling);
        EXPECT_DOUBLE_EQ(H_p(0, 0), 6.0) << "jacobi_scaling = " << jacobi_scaling;
        EXPECT_DOUBLE_EQ(H_p(0, 1), -4.0) << "jacobi_scaling = " << jacobi_scaling;
        EXPECT_DOUBLE_EQ(H_p(1, 0), -4.0) << "jacobi_scaling = " << jacobi_scaling;
        EXPECT_DOUBLE_EQ(H_p(1, 1), 8.0) << "jacobi_scaling = " << jacobi_scaling;
        EXPECT_DOUBLE_EQ(b_p(0), 1.0) << "jacobi_scaling = " << jacobi_scaling;
        EXPECT_DOUBLE_EQ(b_p(1), 2.0) << "jacobi_scaling = " << jacobi_scaling;
    }
}

TEST(schur_complement, jacobi_scaling_equivalence_well_conditioned) {
    // Jacobi scaling is exact, so the paths agree in value. They are not bitwise identical, because the scaling
    // introduces roundings and changes the order of operations.
    const Eigen::MatrixXd H = well_conditioned_spd(6);
    const Eigen::VectorXd b = deterministic_vector(6);
    Eigen::MatrixXd H_p_unscaled, H_p_scaled;
    Eigen::VectorXd b_p_unscaled, b_p_scaled;
    math::schur_complement(H, b, 3, H_p_unscaled, b_p_unscaled, 0.0, 1.0e-9, false);
    math::schur_complement(H, b, 3, H_p_scaled, b_p_scaled, 0.0, 1.0e-9, true);

    EXPECT_TRUE(H_p_scaled.isApprox(H_p_unscaled));
    EXPECT_TRUE(b_p_scaled.isApprox(b_p_unscaled));
}

TEST(schur_complement, jacobi_scaling_equivalence_with_damping) {
    // Scaling the damped block is equivalent to adding damping / H_mm(i, i) to the unit diagonal, so damping keeps
    // its unscaled meaning and the paths agree for a non-zero damping factor.
    const Eigen::MatrixXd H = well_conditioned_spd(6);
    const Eigen::VectorXd b = deterministic_vector(6);
    Eigen::MatrixXd H_p_unscaled, H_p_scaled;
    Eigen::VectorXd b_p_unscaled, b_p_scaled;
    const double damping_unscaled = math::schur_complement(H, b, 3, H_p_unscaled, b_p_unscaled, 0.1, 1.0e-9, false);
    const double damping_scaled = math::schur_complement(H, b, 3, H_p_scaled, b_p_scaled, 0.1, 1.0e-9, true);

    EXPECT_DOUBLE_EQ(damping_scaled, damping_unscaled);
    EXPECT_GT(damping_scaled, 0.0);
    EXPECT_TRUE(H_p_scaled.isApprox(H_p_unscaled));
    EXPECT_TRUE(b_p_scaled.isApprox(b_p_unscaled));
}

TEST(schur_complement, jacobi_scaling_zero_diagonal_throws) {
    // An unconstrained marginalised variable must still fail in the decomposition, rather than the scaling turning it
    // into an infinity or NaN that propagates silently.
    Eigen::MatrixXd H(3, 3);
    H << 0.0, 0.0, 0.0, 0.0, 4.0, 1.0, 0.0, 1.0, 5.0;
    Eigen::VectorXd b(3);
    b << 1.0, 2.0, 3.0;

    for (const bool jacobi_scaling : {false, true}) {
        Eigen::MatrixXd H_p;
        Eigen::VectorXd b_p;
        EXPECT_ANY_THROW(math::schur_complement(H, b, 1, H_p, b_p, 0.0, 1.0e-9, jacobi_scaling))
                << "jacobi_scaling = " << jacobi_scaling;
    }
}

TEST(schur_complement, jacobi_scaling_negative_diagonal_throws) {
    // A numerically indefinite block must also fail in the decomposition, which requires the scaling to leave a
    // non-positive diagonal entry unscaled rather than taking the square root of it.
    Eigen::MatrixXd H(3, 3);
    H << -1.0, 0.0, 2.0, 0.0, 4.0, 1.0, 2.0, 1.0, 5.0;
    Eigen::VectorXd b(3);
    b << 1.0, 2.0, 3.0;

    for (const bool jacobi_scaling : {false, true}) {
        Eigen::MatrixXd H_p;
        Eigen::VectorXd b_p;
        EXPECT_ANY_THROW(math::schur_complement(H, b, 1, H_p, b_p, 0.0, 1.0e-9, jacobi_scaling))
                << "jacobi_scaling = " << jacobi_scaling;
    }
}

TEST(schur_complement, jacobi_scaling_asymmetric_input_still_detected) {
    // Regression test: scaling the off-diagonal blocks independently, rather than transposing one of them, keeps an
    // asymmetry of H detectable in H_p. Transposing would make H_p symmetric by construction.
    Eigen::MatrixXd H(3, 3);
    H << 4.0, 1.0, 2.0, 1.0 + 1.0e-3, 5.0, 1.0, 2.0 + 1.0e-3, 1.0, 6.0;
    Eigen::VectorXd b(3);
    b << 1.0, 2.0, 3.0;

    for (const bool jacobi_scaling : {false, true}) {
        Eigen::MatrixXd H_p;
        Eigen::VectorXd b_p;
        EXPECT_ANY_THROW(math::schur_complement(H, b, 1, H_p, b_p, 0.0, 1.0e-9, jacobi_scaling))
                << "jacobi_scaling = " << jacobi_scaling;
    }
}

// Deterministic full-rank matrix with entries sin((i + 1) (j + 1)), unlike sin(a i + b j + c), whose angle sum
// identity makes any such matrix rank 2
Eigen::MatrixXd deterministic_matrix(const int rows, const int cols) {
    Eigen::MatrixXd M(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            M(i, j) = std::sin(static_cast<double>((i + 1) * (j + 1)));
        }
    }
    return M;
}

// Exact marginal of 0.5 || J dx + r ||^2 over its first left_block_size variables, for any rank of J_m, by projecting
// out the range of J_m: H_p = J_k^T (I - J_m J_m^+) J_k and b_p = -J_k^T (I - J_m J_m^+) r
void projected_marginal(const Eigen::MatrixXd& J, const Eigen::VectorXd& r, const int left_block_size,
        Eigen::MatrixXd& H_p, Eigen::VectorXd& b_p) {
    const Eigen::MatrixXd J_m = J.leftCols(left_block_size);
    const Eigen::MatrixXd J_k = J.rightCols(J.cols() - left_block_size);
    const Eigen::MatrixXd projector =
            Eigen::MatrixXd::Identity(J.rows(), J.rows()) - J_m * J_m.completeOrthogonalDecomposition().pseudoInverse();
    H_p = J_k.transpose() * projector * J_k;
    b_p = -J_k.transpose() * projector * r;
}

TEST(square_root_schur_complement, single_variables_analytic) {
    // H = J^T J = [2, 0; 0, 6] and b = -J^T r = (-3, -5), so marginalising the first variable gives H_p = 6, b_p = -5
    Eigen::MatrixXd J(3, 2);
    J << 1.0, 1.0, 1.0, -1.0, 0.0, 2.0;
    Eigen::VectorXd r(3);
    r << 1.0, 2.0, 3.0;
    Eigen::MatrixXd S_p;
    Eigen::VectorXd r_p, pivots_mm, pivots_p;
    EXPECT_EQ(math::square_root_schur_complement(J, r, 1, S_p, r_p, pivots_mm, pivots_p, 1.0e3), 1);

    ASSERT_EQ(S_p.rows(), 1);
    ASSERT_EQ(S_p.cols(), 1);
    EXPECT_NEAR(S_p(0, 0) * S_p(0, 0), 6.0, 1.0e-12);
    EXPECT_NEAR(-S_p(0, 0) * r_p(0), -5.0, 1.0e-12);
    EXPECT_NEAR(pivots_mm(0), std::sqrt(2.0), 1.0e-12);
    EXPECT_NEAR(pivots_p(0), std::sqrt(6.0), 1.0e-12);
}

TEST(square_root_schur_complement, equivalent_to_schur_complement) {
    const Eigen::MatrixXd J = deterministic_matrix(12, 7);
    const Eigen::VectorXd r = deterministic_vector(12);
    Eigen::MatrixXd S_p, H_p;
    Eigen::VectorXd r_p, pivots_mm, pivots_p, b_p;
    EXPECT_EQ(math::square_root_schur_complement(J, r, 3, S_p, r_p, pivots_mm, pivots_p, 1.0e3), 3);
    math::schur_complement(J.transpose() * J, -J.transpose() * r, 3, H_p, b_p, 0.0, 1.0e-9, false);

    ASSERT_EQ(S_p.rows(), 4);
    ASSERT_EQ(S_p.cols(), 4);
    EXPECT_TRUE((S_p.transpose() * S_p).isApprox(H_p, 1.0e-12));
    EXPECT_TRUE((-S_p.transpose() * r_p).isApprox(b_p, 1.0e-12));
}

TEST(square_root_schur_complement, rank_deficient_left_block_without_damping) {
    // A duplicated column leaves a direction of dx_m unconstrained, so H_mm is singular
    Eigen::MatrixXd J = deterministic_matrix(10, 6);
    J.col(1) = J.col(0);
    const Eigen::VectorXd r = deterministic_vector(10);
    Eigen::MatrixXd S_p, H_p;
    Eigen::VectorXd r_p, pivots_mm, pivots_p, b_p;
    EXPECT_EQ(math::square_root_schur_complement(J, r, 3, S_p, r_p, pivots_mm, pivots_p, 1.0e3), 2);
    projected_marginal(J, r, 3, H_p, b_p);

    EXPECT_TRUE((S_p.transpose() * S_p).isApprox(H_p, 1.0e-12));
    EXPECT_TRUE((-S_p.transpose() * r_p).isApprox(b_p, 1.0e-12));
}

TEST(square_root_schur_complement, fewer_rows_than_right_block) {
    const Eigen::MatrixXd J = deterministic_matrix(5, 6);
    const Eigen::VectorXd r = deterministic_vector(5);
    Eigen::MatrixXd S_p, H_p;
    Eigen::VectorXd r_p, pivots_mm, pivots_p, b_p;
    EXPECT_EQ(math::square_root_schur_complement(J, r, 2, S_p, r_p, pivots_mm, pivots_p, 1.0e3), 2);
    projected_marginal(J, r, 2, H_p, b_p);

    ASSERT_EQ(S_p.rows(), 3);
    ASSERT_EQ(S_p.cols(), 4);
    ASSERT_EQ(pivots_p.size(), 3);
    EXPECT_TRUE((S_p.transpose() * S_p).isApprox(H_p, 1.0e-12));
    EXPECT_TRUE((-S_p.transpose() * r_p).isApprox(b_p, 1.0e-12));
}

TEST(square_root_schur_complement, all_rows_absorbed_by_left_block) {
    const Eigen::MatrixXd J = deterministic_matrix(2, 3);
    const Eigen::VectorXd r = deterministic_vector(2);
    Eigen::MatrixXd S_p;
    Eigen::VectorXd r_p, pivots_mm, pivots_p;
    EXPECT_EQ(math::square_root_schur_complement(J, r, 2, S_p, r_p, pivots_mm, pivots_p, 1.0e3), 2);

    EXPECT_EQ(S_p.rows(), 0);
    EXPECT_EQ(S_p.cols(), 1);
    EXPECT_EQ(r_p.size(), 0);
    EXPECT_EQ(pivots_p.size(), 0);
}

TEST(square_root_schur_complement, nullspace_rows_removed_by_pivot_threshold) {
    // A duplicated kept column gives the prior a nullspace direction, whose row has a numerically zero pivot
    Eigen::MatrixXd J = deterministic_matrix(12, 7);
    J.col(6) = J.col(5);
    const Eigen::VectorXd r = deterministic_vector(12);
    Eigen::MatrixXd S_p;
    Eigen::VectorXd r_p, pivots_mm, pivots_p;
    math::square_root_schur_complement(J, r, 3, S_p, r_p, pivots_mm, pivots_p, 1.0e3);
    const double nullspace_threshold = 1.0e3 * std::numeric_limits<double>::epsilon() * pivots_p.maxCoeff();
    const Eigen::MatrixXd S_p_reduced = math::remove_rows_by_threshold(S_p, pivots_p, nullspace_threshold);

    ASSERT_EQ(S_p.rows(), 4);
    EXPECT_EQ(S_p_reduced.rows(), 3);
    EXPECT_TRUE((S_p_reduced.transpose() * S_p_reduced).isApprox(S_p.transpose() * S_p, 1.0e-12));
}

TEST(square_root_schur_complement, preserves_weak_information_lost_by_normal_equations) {
    // A stiff relative residual A (x_k - x_m) and a weak absolute residual w x_m, like positions pinned tightly to each
    // other (e.g. by an accelerometer through a trajectory) but only weakly to the world (a gauge direction). The exact
    // marginal information of x_k is A^2 w^2 / (A^2 + w^2), but forming J^T J rounds A^2 + w^2 to A^2, so that the
    // Schur complement of the normal equations loses it entirely.
    const double A = 1.0e7;
    const double w = std::sqrt(1.0e-3);
    Eigen::MatrixXd J(2, 2);
    J << -A, A, w, 0.0;
    const Eigen::VectorXd r = Eigen::VectorXd::Zero(2);
    const double information = A * A * w * w / (A * A + w * w);
    Eigen::MatrixXd S_p, H_p;
    Eigen::VectorXd r_p, pivots_mm, pivots_p, b_p;
    math::square_root_schur_complement(J, r, 1, S_p, r_p, pivots_mm, pivots_p, 1.0e3);
    math::schur_complement(J.transpose() * J, -J.transpose() * r, 1, H_p, b_p, 0.0, 1.0e-9, false);

    EXPECT_NEAR(S_p(0, 0) * S_p(0, 0), information, 1.0e-12 * information);
    EXPECT_LT(H_p(0, 0), 0.5 * information);
}

TEST(reorder_symmetric_matrix, index_1_3x3) {
    Eigen::Matrix3d cov;
    cov << 1.0, 4.0, 5.0, 4.0, 2.0, 6.0, 5.0, 6.0, 3.0;
    Eigen::Matrix3d reordered_cov_truth;
    reordered_cov_truth << 2.0, 6.0, 4.0, 6.0, 3.0, 5.0, 4.0, 5.0, 1.0;
    const Eigen::Matrix3d reordered_cov = math::reorder_symmetric_matrix(cov, 1);
    EXPECT_TRUE(reordered_cov.isApprox(reordered_cov_truth));
}

TEST(reorder_symmetric_matrix, index_1_6x6) {
    Eigen::Matrix<double, 6, 6> cov;
    cov << 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.002, 0.007, 0.008, 0.009, 0.010, 0.011, 0.003, 0.008, 0.012,
            0.013, 0.014, 0.015, 0.004, 0.009, 0.013, 0.016, 0.017, 0.018, 0.005, 0.010, 0.014, 0.017, 0.019, 0.020,
            0.006, 0.011, 0.015, 0.018, 0.020, 0.021;
    Eigen::Matrix<double, 6, 6> reordered_cov_truth;
    reordered_cov_truth << 0.007, 0.008, 0.009, 0.010, 0.011, 0.002, 0.008, 0.012, 0.013, 0.014, 0.015, 0.003, 0.009,
            0.013, 0.016, 0.017, 0.018, 0.004, 0.010, 0.014, 0.017, 0.019, 0.020, 0.005, 0.011, 0.015, 0.018, 0.020,
            0.021, 0.006, 0.002, 0.003, 0.004, 0.005, 0.006, 0.001;
    Eigen::Matrix<double, 6, 6> reordered_cov = math::reorder_symmetric_matrix(cov, 1);
    EXPECT_TRUE(reordered_cov.isApprox(reordered_cov_truth));
}

TEST(reorder_symmetric_matrix, index_2_6x6) {
    Eigen::Matrix<double, 6, 6> cov;
    cov << 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.002, 0.007, 0.008, 0.009, 0.010, 0.011, 0.003, 0.008, 0.012,
            0.013, 0.014, 0.015, 0.004, 0.009, 0.013, 0.016, 0.017, 0.018, 0.005, 0.010, 0.014, 0.017, 0.019, 0.020,
            0.006, 0.011, 0.015, 0.018, 0.020, 0.021;
    Eigen::Matrix<double, 6, 6> reordered_cov_truth;
    reordered_cov_truth << 0.012, 0.013, 0.014, 0.015, 0.003, 0.008, 0.013, 0.016, 0.017, 0.018, 0.004, 0.009, 0.014,
            0.017, 0.019, 0.020, 0.005, 0.010, 0.015, 0.018, 0.020, 0.021, 0.006, 0.011, 0.003, 0.004, 0.005, 0.006,
            0.001, 0.002, 0.008, 0.009, 0.010, 0.011, 0.002, 0.007;
    Eigen::Matrix<double, 6, 6> reordered_cov = math::reorder_symmetric_matrix(cov, 2);
    EXPECT_TRUE(reordered_cov.isApprox(reordered_cov_truth));
}

TEST(reorder_symmetric_matrix, index_3_6x6) {
    Eigen::Matrix<double, 6, 6> cov;
    cov << 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.002, 0.007, 0.008, 0.009, 0.010, 0.011, 0.003, 0.008, 0.012,
            0.013, 0.014, 0.015, 0.004, 0.009, 0.013, 0.016, 0.017, 0.018, 0.005, 0.010, 0.014, 0.017, 0.019, 0.020,
            0.006, 0.011, 0.015, 0.018, 0.020, 0.021;
    Eigen::Matrix<double, 6, 6> reordered_cov_truth;
    reordered_cov_truth << 0.016, 0.017, 0.018, 0.004, 0.009, 0.013, 0.017, 0.019, 0.020, 0.005, 0.010, 0.014, 0.018,
            0.020, 0.021, 0.006, 0.011, 0.015, 0.004, 0.005, 0.006, 0.001, 0.002, 0.003, 0.009, 0.010, 0.011, 0.002,
            0.007, 0.008, 0.013, 0.014, 0.015, 0.003, 0.008, 0.012;
    Eigen::Matrix<double, 6, 6> reordered_cov = math::reorder_symmetric_matrix(cov, 3);
    EXPECT_TRUE(reordered_cov.isApprox(reordered_cov_truth));
}

TEST(reorder_symmetric_matrix, index_4_6x6) {
    Eigen::Matrix<double, 6, 6> cov;
    cov << 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.002, 0.007, 0.008, 0.009, 0.010, 0.011, 0.003, 0.008, 0.012,
            0.013, 0.014, 0.015, 0.004, 0.009, 0.013, 0.016, 0.017, 0.018, 0.005, 0.010, 0.014, 0.017, 0.019, 0.020,
            0.006, 0.011, 0.015, 0.018, 0.020, 0.021;
    Eigen::Matrix<double, 6, 6> reordered_cov_truth;
    reordered_cov_truth << 0.019, 0.020, 0.005, 0.010, 0.014, 0.017, 0.020, 0.021, 0.006, 0.011, 0.015, 0.018, 0.005,
            0.006, 0.001, 0.002, 0.003, 0.004, 0.010, 0.011, 0.002, 0.007, 0.008, 0.009, 0.014, 0.015, 0.003, 0.008,
            0.012, 0.013, 0.017, 0.018, 0.004, 0.009, 0.013, 0.016;
    Eigen::Matrix<double, 6, 6> reordered_cov = math::reorder_symmetric_matrix(cov, 4);
    EXPECT_TRUE(reordered_cov.isApprox(reordered_cov_truth));
}

TEST(reorder_symmetric_matrix, index_5_6x6) {
    Eigen::Matrix<double, 6, 6> cov;
    cov << 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.002, 0.007, 0.008, 0.009, 0.010, 0.011, 0.003, 0.008, 0.012,
            0.013, 0.014, 0.015, 0.004, 0.009, 0.013, 0.016, 0.017, 0.018, 0.005, 0.010, 0.014, 0.017, 0.019, 0.020,
            0.006, 0.011, 0.015, 0.018, 0.020, 0.021;
    Eigen::Matrix<double, 6, 6> reordered_cov_truth;
    reordered_cov_truth << 0.021, 0.006, 0.011, 0.015, 0.018, 0.020, 0.006, 0.001, 0.002, 0.003, 0.004, 0.005, 0.011,
            0.002, 0.007, 0.008, 0.009, 0.010, 0.015, 0.003, 0.008, 0.012, 0.013, 0.014, 0.018, 0.004, 0.009, 0.013,
            0.016, 0.017, 0.020, 0.005, 0.010, 0.014, 0.017, 0.019;
    Eigen::Matrix<double, 6, 6> reordered_cov = math::reorder_symmetric_matrix(cov, 5);
    EXPECT_TRUE(reordered_cov.isApprox(reordered_cov_truth));
}
