#include <gtest/gtest.h>

#include <Eigen/Core>

#include "mathbox/stiffness.hpp"

TEST(decompose, LLT_cholesky_I) {
    const Eigen::Matrix3d matrix = Eigen::Matrix3d::Identity();
    const Eigen::Matrix3d L = math::LLT(matrix, math::LLTDecompositionMethod::CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(L * L.transpose()));
}

TEST(decompose, LLT_cholesky_ascending) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 1.0, 2.0, 3.0).finished().asDiagonal();
    const Eigen::Matrix3d L = math::LLT(matrix, math::LLTDecompositionMethod::CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(L * L.transpose()));
}

TEST(decompose, LLT_cholesky_descending) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 3.0, 2.0, 1.0).finished().asDiagonal();
    const Eigen::Matrix3d L = math::LLT(matrix, math::LLTDecompositionMethod::CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(L * L.transpose()));
}

TEST(decompose, LLT_cholesky_mixed_1) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 1.0, 3.0, 2.0).finished().asDiagonal();
    const Eigen::Matrix3d L = math::LLT(matrix, math::LLTDecompositionMethod::CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(L * L.transpose()));
}

TEST(decompose, LLT_cholesky_mixed_2) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 2.0, 1.0, 3.0).finished().asDiagonal();
    const Eigen::Matrix3d L = math::LLT(matrix, math::LLTDecompositionMethod::CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(L * L.transpose()));
}

TEST(decompose, LLT_cholesky_mixed_3) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 2.0, 3.0, 1.0).finished().asDiagonal();
    const Eigen::Matrix3d L = math::LLT(matrix, math::LLTDecompositionMethod::CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(L * L.transpose()));
}

TEST(decompose, LLT_cholesky_mixed_4) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 3.0, 1.0, 2.0).finished().asDiagonal();
    const Eigen::Matrix3d L = math::LLT(matrix, math::LLTDecompositionMethod::CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(L * L.transpose()));
}

TEST(decompose, LLT_cholesky_symmetric_positive_definite) {
    const Eigen::Matrix3d matrix = (Eigen::Matrix3d() << 2.5, 0.4, -0.9, 0.4, 0.8, 0.9, -0.9, 0.9, 6.3).finished();
    const Eigen::Matrix3d L = math::LLT(matrix, math::LLTDecompositionMethod::CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(L * L.transpose()));
}

TEST(decompose, LLT_eigen_I) {
    const Eigen::Matrix3d matrix = Eigen::Matrix3d::Identity();
    const Eigen::Matrix3d L = math::LLT(matrix, math::LLTDecompositionMethod::EIGEN);
    EXPECT_TRUE(matrix.isApprox(L * L.transpose()));
}

TEST(decompose, LLT_eigen_ascending) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 1.0, 2.0, 3.0).finished().asDiagonal();
    const Eigen::Matrix3d L = math::LLT(matrix, math::LLTDecompositionMethod::EIGEN);
    EXPECT_TRUE(matrix.isApprox(L * L.transpose()));
}

TEST(decompose, LLT_eigen_descending) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 3.0, 2.0, 1.0).finished().asDiagonal();
    const Eigen::Matrix3d L = math::LLT(matrix, math::LLTDecompositionMethod::EIGEN);
    EXPECT_TRUE(matrix.isApprox(L * L.transpose()));
}

TEST(decompose, LLT_eigen_mixed_1) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 1.0, 3.0, 2.0).finished().asDiagonal();
    const Eigen::Matrix3d L = math::LLT(matrix, math::LLTDecompositionMethod::EIGEN);
    EXPECT_TRUE(matrix.isApprox(L * L.transpose()));
}

TEST(decompose, LLT_eigen_mixed_2) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 2.0, 1.0, 3.0).finished().asDiagonal();
    const Eigen::Matrix3d L = math::LLT(matrix, math::LLTDecompositionMethod::EIGEN);
    EXPECT_TRUE(matrix.isApprox(L * L.transpose()));
}

TEST(decompose, LLT_eigen_mixed_3) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 2.0, 3.0, 1.0).finished().asDiagonal();
    const Eigen::Matrix3d L = math::LLT(matrix, math::LLTDecompositionMethod::EIGEN);
    EXPECT_TRUE(matrix.isApprox(L * L.transpose()));
}

TEST(decompose, LLT_eigen_mixed_4) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 3.0, 1.0, 2.0).finished().asDiagonal();
    const Eigen::Matrix3d L = math::LLT(matrix, math::LLTDecompositionMethod::EIGEN);
    EXPECT_TRUE(matrix.isApprox(L * L.transpose()));
}

TEST(decompose, LLT_eigen_symmetric_positive_definite) {
    const Eigen::Matrix3d matrix = (Eigen::Matrix3d() << 2.5, 0.4, -0.9, 0.4, 0.8, 0.9, -0.9, 0.9, 6.3).finished();
    const Eigen::Matrix3d L = math::LLT(matrix, math::LLTDecompositionMethod::EIGEN);
    EXPECT_TRUE(matrix.isApprox(L * L.transpose()));
}

TEST(decompose, LLT_robust_cholesky_I) {
    const Eigen::Matrix3d matrix = Eigen::Matrix3d::Identity();
    const Eigen::Matrix3d L = math::LLT(matrix, math::LLTDecompositionMethod::ROBUST_CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(L * L.transpose()));
}

TEST(decompose, LLT_robust_cholesky_ascending) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 1.0, 2.0, 3.0).finished().asDiagonal();
    const Eigen::Matrix3d L = math::LLT(matrix, math::LLTDecompositionMethod::ROBUST_CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(L * L.transpose()));
}

TEST(decompose, LLT_robust_cholesky_descending) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 3.0, 2.0, 1.0).finished().asDiagonal();
    const Eigen::Matrix3d L = math::LLT(matrix, math::LLTDecompositionMethod::ROBUST_CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(L * L.transpose()));
}

TEST(decompose, LLT_robust_cholesky_mixed_1) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 1.0, 3.0, 2.0).finished().asDiagonal();
    const Eigen::Matrix3d L = math::LLT(matrix, math::LLTDecompositionMethod::ROBUST_CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(L * L.transpose()));
}

TEST(decompose, LLT_robust_cholesky_mixed_2) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 2.0, 1.0, 3.0).finished().asDiagonal();
    const Eigen::Matrix3d L = math::LLT(matrix, math::LLTDecompositionMethod::ROBUST_CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(L * L.transpose()));
}

TEST(decompose, LLT_robust_cholesky_mixed_3) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 2.0, 3.0, 1.0).finished().asDiagonal();
    const Eigen::Matrix3d L = math::LLT(matrix, math::LLTDecompositionMethod::ROBUST_CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(L * L.transpose()));
}

TEST(decompose, LLT_robust_cholesky_mixed_4) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 3.0, 1.0, 2.0).finished().asDiagonal();
    const Eigen::Matrix3d L = math::LLT(matrix, math::LLTDecompositionMethod::ROBUST_CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(L * L.transpose()));
}

TEST(decompose, LLT_robust_cholesky_symmetric_positive_definite) {
    const Eigen::Matrix3d matrix = (Eigen::Matrix3d() << 2.5, 0.4, -0.9, 0.4, 0.8, 0.9, -0.9, 0.9, 6.3).finished();
    const Eigen::Matrix3d L = math::LLT(matrix, math::LLTDecompositionMethod::ROBUST_CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(L * L.transpose()));
}
TEST(decompose, UTU_cholesky_I) {
    const Eigen::Matrix3d matrix = Eigen::Matrix3d::Identity();
    const Eigen::Matrix3d U = math::UTU(matrix, math::LLTDecompositionMethod::CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(U.transpose() * U));
}

TEST(decompose, UTU_cholesky_ascending) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 1.0, 2.0, 3.0).finished().asDiagonal();
    const Eigen::Matrix3d U = math::UTU(matrix, math::LLTDecompositionMethod::CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(U.transpose() * U));
}

TEST(decompose, UTU_cholesky_descending) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 3.0, 2.0, 1.0).finished().asDiagonal();
    const Eigen::Matrix3d U = math::UTU(matrix, math::LLTDecompositionMethod::CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(U.transpose() * U));
}

TEST(decompose, UTU_cholesky_mixed_1) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 1.0, 3.0, 2.0).finished().asDiagonal();
    const Eigen::Matrix3d U = math::UTU(matrix, math::LLTDecompositionMethod::CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(U.transpose() * U));
}

TEST(decompose, UTU_cholesky_mixed_2) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 2.0, 1.0, 3.0).finished().asDiagonal();
    const Eigen::Matrix3d U = math::UTU(matrix, math::LLTDecompositionMethod::CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(U.transpose() * U));
}

TEST(decompose, UTU_cholesky_mixed_3) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 2.0, 3.0, 1.0).finished().asDiagonal();
    const Eigen::Matrix3d U = math::UTU(matrix, math::LLTDecompositionMethod::CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(U.transpose() * U));
}

TEST(decompose, UTU_cholesky_mixed_4) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 3.0, 1.0, 2.0).finished().asDiagonal();
    const Eigen::Matrix3d U = math::UTU(matrix, math::LLTDecompositionMethod::CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(U.transpose() * U));
}

TEST(decompose, UTU_cholesky_symmetric_positive_definite) {
    const Eigen::Matrix3d matrix = (Eigen::Matrix3d() << 2.5, 0.4, -0.9, 0.4, 0.8, 0.9, -0.9, 0.9, 6.3).finished();
    const Eigen::Matrix3d U = math::UTU(matrix, math::LLTDecompositionMethod::CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(U.transpose() * U));
}

TEST(decompose, UTU_eigen_I) {
    const Eigen::Matrix3d matrix = Eigen::Matrix3d::Identity();
    const Eigen::Matrix3d U = math::UTU(matrix, math::LLTDecompositionMethod::EIGEN);
    EXPECT_TRUE(matrix.isApprox(U.transpose() * U));
}

TEST(decompose, UTU_eigen_ascending) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 1.0, 2.0, 3.0).finished().asDiagonal();
    const Eigen::Matrix3d U = math::UTU(matrix, math::LLTDecompositionMethod::EIGEN);
    EXPECT_TRUE(matrix.isApprox(U.transpose() * U));
}

TEST(decompose, UTU_eigen_descending) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 3.0, 2.0, 1.0).finished().asDiagonal();
    const Eigen::Matrix3d U = math::UTU(matrix, math::LLTDecompositionMethod::EIGEN);
    EXPECT_TRUE(matrix.isApprox(U.transpose() * U));
}

TEST(decompose, UTU_eigen_mixed_1) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 1.0, 3.0, 2.0).finished().asDiagonal();
    const Eigen::Matrix3d U = math::UTU(matrix, math::LLTDecompositionMethod::EIGEN);
    EXPECT_TRUE(matrix.isApprox(U.transpose() * U));
}

TEST(decompose, UTU_eigen_mixed_2) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 2.0, 1.0, 3.0).finished().asDiagonal();
    const Eigen::Matrix3d U = math::UTU(matrix, math::LLTDecompositionMethod::EIGEN);
    EXPECT_TRUE(matrix.isApprox(U.transpose() * U));
}

TEST(decompose, UTU_eigen_mixed_3) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 2.0, 3.0, 1.0).finished().asDiagonal();
    const Eigen::Matrix3d U = math::UTU(matrix, math::LLTDecompositionMethod::EIGEN);
    EXPECT_TRUE(matrix.isApprox(U.transpose() * U));
}

TEST(decompose, UTU_eigen_mixed_4) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 3.0, 1.0, 2.0).finished().asDiagonal();
    const Eigen::Matrix3d U = math::UTU(matrix, math::LLTDecompositionMethod::EIGEN);
    EXPECT_TRUE(matrix.isApprox(U.transpose() * U));
}

TEST(decompose, UTU_eigen_symmetric_positive_definite) {
    const Eigen::Matrix3d matrix = (Eigen::Matrix3d() << 2.5, 0.4, -0.9, 0.4, 0.8, 0.9, -0.9, 0.9, 6.3).finished();
    const Eigen::Matrix3d U = math::UTU(matrix, math::LLTDecompositionMethod::EIGEN);
    EXPECT_TRUE(matrix.isApprox(U.transpose() * U));
}

TEST(decompose, UTU_robust_cholesky_I) {
    const Eigen::Matrix3d matrix = Eigen::Matrix3d::Identity();
    const Eigen::Matrix3d U = math::UTU(matrix, math::LLTDecompositionMethod::ROBUST_CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(U.transpose() * U));
}

TEST(decompose, UTU_robust_cholesky_ascending) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 1.0, 2.0, 3.0).finished().asDiagonal();
    const Eigen::Matrix3d U = math::UTU(matrix, math::LLTDecompositionMethod::ROBUST_CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(U.transpose() * U));
}

TEST(decompose, UTU_robust_cholesky_descending) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 3.0, 2.0, 1.0).finished().asDiagonal();
    const Eigen::Matrix3d U = math::UTU(matrix, math::LLTDecompositionMethod::ROBUST_CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(U.transpose() * U));
}

TEST(decompose, UTU_robust_cholesky_mixed_1) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 1.0, 3.0, 2.0).finished().asDiagonal();
    const Eigen::Matrix3d U = math::UTU(matrix, math::LLTDecompositionMethod::ROBUST_CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(U.transpose() * U));
}

TEST(decompose, UTU_robust_cholesky_mixed_2) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 2.0, 1.0, 3.0).finished().asDiagonal();
    const Eigen::Matrix3d U = math::UTU(matrix, math::LLTDecompositionMethod::ROBUST_CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(U.transpose() * U));
}

TEST(decompose, UTU_robust_cholesky_mixed_3) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 2.0, 3.0, 1.0).finished().asDiagonal();
    const Eigen::Matrix3d U = math::UTU(matrix, math::LLTDecompositionMethod::ROBUST_CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(U.transpose() * U));
}

TEST(decompose, UTU_robust_cholesky_mixed_4) {
    const Eigen::Matrix3d matrix = (Eigen::Vector3d() << 3.0, 1.0, 2.0).finished().asDiagonal();
    const Eigen::Matrix3d U = math::UTU(matrix, math::LLTDecompositionMethod::ROBUST_CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(U.transpose() * U));
}

TEST(decompose, UTU_robust_cholesky_symmetric_positive_definite) {
    const Eigen::Matrix3d matrix = (Eigen::Matrix3d() << 2.5, 0.4, -0.9, 0.4, 0.8, 0.9, -0.9, 0.9, 6.3).finished();
    const Eigen::Matrix3d U = math::UTU(matrix, math::LLTDecompositionMethod::ROBUST_CHOLESKY);
    EXPECT_TRUE(matrix.isApprox(U.transpose() * U));
}
