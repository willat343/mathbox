#include "mathbox/matrix_operations.hpp"

#include <gtest/gtest.h>

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
