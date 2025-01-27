#include "mathbox/covariance.hpp"

#include <gtest/gtest.h>

TEST(covariance, constructors_fixed) {
    math::Covariance3d covariance0(0.1 * Eigen::Matrix3d::Identity());
    math::Covariance3d covariance1((Eigen::Matrix3d() << 0.1, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.1).finished());
    math::Covariance3d covariance2((Eigen::Vector3d() << 0.1, 0.1, 0.1).finished());
    math::Covariance3d covariance3(0.1, 3);
    math::Covariance3d covariance4(0.1);
    EXPECT_EQ(covariance0, covariance1);
    EXPECT_EQ(covariance0, covariance2);
    EXPECT_EQ(covariance0, covariance3);
    EXPECT_EQ(covariance0, covariance4);
}

TEST(covariance, constructors_dynamic) {
    math::CovarianceXd covariance0(0.1 * Eigen::MatrixXd::Identity(3, 3));
    math::CovarianceXd covariance1((Eigen::MatrixXd(3, 3) << 0.1, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.1).finished());
    math::CovarianceXd covariance2((Eigen::VectorXd(3) << 0.1, 0.1, 0.1).finished());
    math::CovarianceXd covariance3(0.1, 3);
    EXPECT_EQ(covariance0, covariance1);
    EXPECT_EQ(covariance0, covariance2);
    EXPECT_EQ(covariance0, covariance3);
}

TEST(covariance_density, covariance_fixed) {
    math::CovarianceDensity3d covariance_density(0.1 * Eigen::Matrix3d::Identity());
    const double span = 3.4;
    math::Covariance3d covariance(span * 0.1 * Eigen::Matrix3d::Identity());
    EXPECT_EQ(covariance_density.covariance(span), covariance);
}

TEST(covariance_density, covariance_dynamic) {
    math::CovarianceDensityXd covariance_density(0.1 * Eigen::MatrixXd::Identity(3, 3));
    const double span = 3.4;
    math::CovarianceXd covariance(span * 0.1 * Eigen::MatrixXd::Identity(3, 3));
    EXPECT_EQ(covariance_density.covariance(span), covariance);
}
