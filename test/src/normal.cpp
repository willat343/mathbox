#include "mathbox/normal.hpp"

#include <gtest/gtest.h>

#include <optional>

TEST(GrvGenerator, fixed_and_dynamic_1D) {
    for (const auto method : math::llt_decomposition_methods) {
        math::GrvGenerator<double, 1> grvg_fixed(Eigen::Matrix<double, 1, 1>(3.0),
                (Eigen::Matrix<double, 1, 1>() << 2.3).finished(), method, 0);
        math::GrvGenerator<double, Eigen::Dynamic> grvg_dynamic((Eigen::VectorXd(1) << 3.0).finished(),
                (Eigen::MatrixXd(1, 1) << 2.3).finished(), method, 0);
        EXPECT_TRUE(grvg_fixed.compute_covariance().isApprox(grvg_dynamic.compute_covariance()));
        EXPECT_TRUE(grvg_fixed().isApprox(grvg_dynamic()));
    }
}

TEST(GrvGenerator, fixed_and_dynamic_2D) {
    for (const auto method : math::llt_decomposition_methods) {
        math::GrvGenerator<double, 2> grvg_fixed(Eigen::Vector2d(2.0, 3.0),
                (Eigen::Matrix2d() << 2.3, 0.3, 0.3, 1.2).finished(), method, 0);
        math::GrvGenerator<double, Eigen::Dynamic> grvg_dynamic((Eigen::VectorXd(2) << 2.0, 3.0).finished(),
                (Eigen::MatrixXd(2, 2) << 2.3, 0.3, 0.3, 1.2).finished(), method, 0);
        EXPECT_TRUE(grvg_fixed.compute_covariance().isApprox(grvg_dynamic.compute_covariance()));
        EXPECT_TRUE(grvg_fixed().isApprox(grvg_dynamic()));
    }
}

TEST(GrvGenerator, fixed_and_dynamic_3D) {
    for (const auto method : math::llt_decomposition_methods) {
        math::GrvGenerator<double, 3> grvg_fixed(Eigen::Vector3d(1.0, 2.0, 3.0),
                (Eigen::Matrix3d() << 2.3, 0.3, -0.4, 0.3, 1.2, 0.2, -0.4, 0.2, 3.3).finished(), method, 0);
        math::GrvGenerator<double, Eigen::Dynamic> grvg_dynamic((Eigen::VectorXd(3) << 1.0, 2.0, 3.0).finished(),
                (Eigen::MatrixXd(3, 3) << 2.3, 0.3, -0.4, 0.3, 1.2, 0.2, -0.4, 0.2, 3.3).finished(), method, 0);
        EXPECT_TRUE(grvg_fixed.compute_covariance().isApprox(grvg_dynamic.compute_covariance()));
        EXPECT_TRUE(grvg_fixed().isApprox(grvg_dynamic()));
    }
}

TEST(GrvGenerator, decomposers_1D) {
    const Eigen::Matrix<double, 1, 1> mean = Eigen::Matrix<double, 1, 1>(3.0);
    const Eigen::Matrix<double, 1, 1> covariance = (Eigen::Matrix<double, 1, 1>() << 2.0).finished();
    std::optional<Eigen::Matrix<double, 1, 1>> transform;
    std::optional<Eigen::Matrix<double, 1, 1>> sample;
    for (const auto method : math::llt_decomposition_methods) {
        math::GrvGenerator<double, 1> grvg(mean, covariance, method, 0);
        EXPECT_TRUE(grvg.compute_covariance().isApprox(covariance));
        // Transform is guaranteed to be the same because there is no ordering ambiguity.
        if (transform) {
            EXPECT_TRUE(grvg.transform().isApprox(transform.value()));
        } else {
            transform = grvg.transform();
        }
        // Sample is guaranteed to be the same because transform is the same.
        if (sample) {
            EXPECT_TRUE(grvg().isApprox(sample.value()));
        } else {
            sample = grvg();
        }
    }
}

TEST(GrvGenerator, decomposers_2D_uncorrelated) {
    const Eigen::Vector2d mean = Eigen::Vector2d(2.0, 3.0);
    const Eigen::Matrix2d covariance = (Eigen::Matrix2d() << 2.0, 0.0, 0.0, 1.5).finished();
    for (const auto method : math::llt_decomposition_methods) {
        math::GrvGenerator<double, 2> grvg(mean, covariance, method, 0);
        EXPECT_TRUE(grvg.compute_covariance().isApprox(covariance));
    }
}

TEST(GrvGenerator, decomposers_2D_correlated) {
    const Eigen::Vector2d mean = Eigen::Vector2d(2.0, 3.0);
    const Eigen::Matrix2d covariance = (Eigen::Matrix2d() << 2.0, 0.3, 0.3, 1.5).finished();
    for (const auto method : math::llt_decomposition_methods) {
        math::GrvGenerator<double, 2> grvg(mean, covariance, method, 0);
        EXPECT_TRUE(grvg.compute_covariance().isApprox(covariance));
    }
}

TEST(GrvGenerator, decomposers_3D_uncorrelated) {
    const Eigen::Vector3d mean = Eigen::Vector3d(1.0, 2.0, 3.0);
    const Eigen::Matrix3d covariance = (Eigen::Matrix3d() << 2.3, 0.0, 0.0, 0.0, 1.2, 0.0, 0.0, 0.0, 3.3).finished();
    for (const auto method : math::llt_decomposition_methods) {
        math::GrvGenerator<double, 3> grvg(mean, covariance, method, 0);
        EXPECT_TRUE(grvg.compute_covariance().isApprox(covariance));
    }
}

TEST(GrvGenerator, decomposers_3D_correlated) {
    const Eigen::Vector3d mean = Eigen::Vector3d(1.0, 2.0, 3.0);
    const Eigen::Matrix3d covariance = (Eigen::Matrix3d() << 2.3, 0.3, -0.4, 0.3, 1.2, 0.2, -0.4, 0.2, 3.3).finished();
    for (const auto method : math::llt_decomposition_methods) {
        math::GrvGenerator<double, 3> grvg(mean, covariance, method, 0);
        EXPECT_TRUE(grvg.compute_covariance().isApprox(covariance));
    }
}
