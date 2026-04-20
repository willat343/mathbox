#include "mathbox/traits.hpp"

#include <gtest/gtest.h>

TEST(traits, ref) {
    EXPECT_TRUE((std::is_same_v<cppbox::ref_t<Eigen::Vector3d>, Eigen::Ref<Eigen::Vector3d>>));
    EXPECT_TRUE((std::is_same_v<cppbox::ref_t<Eigen::VectorXd>, Eigen::Ref<Eigen::VectorXd>>));
    EXPECT_TRUE((std::is_same_v<cppbox::ref_t<Eigen::Matrix3d>, Eigen::Ref<Eigen::Matrix3d>>));
    EXPECT_TRUE((std::is_same_v<cppbox::ref_t<Eigen::MatrixXd>, Eigen::Ref<Eigen::MatrixXd>>));
}

TEST(traits, const_ref) {
    EXPECT_TRUE((std::is_same_v<cppbox::const_ref_t<Eigen::Vector3d>, const Eigen::Ref<const Eigen::Vector3d>&>));
    EXPECT_TRUE((std::is_same_v<cppbox::const_ref_t<Eigen::VectorXd>, const Eigen::Ref<const Eigen::VectorXd>&>));
    EXPECT_TRUE((std::is_same_v<cppbox::const_ref_t<Eigen::Matrix3d>, const Eigen::Ref<const Eigen::Matrix3d>&>));
    EXPECT_TRUE((std::is_same_v<cppbox::const_ref_t<Eigen::MatrixXd>, const Eigen::Ref<const Eigen::MatrixXd>&>));
}

TEST(traits, remove_ref) {
    EXPECT_TRUE((std::is_same_v<cppbox::remove_ref_t<Eigen::Vector3d>, Eigen::Vector3d>));
    EXPECT_TRUE((std::is_same_v<cppbox::remove_ref_t<Eigen::VectorXd>, Eigen::VectorXd>));
    EXPECT_TRUE((std::is_same_v<cppbox::remove_ref_t<Eigen::Matrix3d>, Eigen::Matrix3d>));
    EXPECT_TRUE((std::is_same_v<cppbox::remove_ref_t<Eigen::MatrixXd>, Eigen::MatrixXd>));
    EXPECT_TRUE((std::is_same_v<cppbox::remove_ref_t<Eigen::Vector3d&>, Eigen::Vector3d>));
    EXPECT_TRUE((std::is_same_v<cppbox::remove_ref_t<Eigen::VectorXd&>, Eigen::VectorXd>));
    EXPECT_TRUE((std::is_same_v<cppbox::remove_ref_t<Eigen::Matrix3d&>, Eigen::Matrix3d>));
    EXPECT_TRUE((std::is_same_v<cppbox::remove_ref_t<Eigen::MatrixXd&>, Eigen::MatrixXd>));
    EXPECT_TRUE((std::is_same_v<cppbox::remove_ref_t<Eigen::Ref<Eigen::Vector3d>>, Eigen::Vector3d>));
    EXPECT_TRUE((std::is_same_v<cppbox::remove_ref_t<Eigen::Ref<Eigen::VectorXd>>, Eigen::VectorXd>));
    EXPECT_TRUE((std::is_same_v<cppbox::remove_ref_t<Eigen::Ref<Eigen::Matrix3d>>, Eigen::Matrix3d>));
    EXPECT_TRUE((std::is_same_v<cppbox::remove_ref_t<Eigen::Ref<Eigen::MatrixXd>>, Eigen::MatrixXd>));
    EXPECT_TRUE((std::is_same_v<cppbox::remove_ref_t<Eigen::Ref<const Eigen::Vector3d>>, const Eigen::Vector3d>));
    EXPECT_TRUE((std::is_same_v<cppbox::remove_ref_t<Eigen::Ref<const Eigen::VectorXd>>, const Eigen::VectorXd>));
    EXPECT_TRUE((std::is_same_v<cppbox::remove_ref_t<Eigen::Ref<const Eigen::Matrix3d>>, const Eigen::Matrix3d>));
    EXPECT_TRUE((std::is_same_v<cppbox::remove_ref_t<Eigen::Ref<const Eigen::MatrixXd>>, const Eigen::MatrixXd>));
}

TEST(traits, is_2d_or_3d) {
    EXPECT_FALSE(math::is_2d_or_3d<-1>);
    EXPECT_FALSE(math::is_2d_or_3d<0>);
    EXPECT_FALSE(math::is_2d_or_3d<1>);
    EXPECT_TRUE(math::is_2d_or_3d<2>);
    EXPECT_TRUE(math::is_2d_or_3d<3>);
    EXPECT_FALSE(math::is_2d_or_3d<4>);
    EXPECT_FALSE(math::is_2d_or_3d<5>);
}

TEST(traits, IsMatrix) {
    static_assert(!math::IsMatrix<double>);
    static_assert(!math::IsMatrix<double&>);
    static_assert(math::IsMatrix<Eigen::VectorXd>);
    static_assert(math::IsMatrix<Eigen::Vector<double, 1>>);
    static_assert(math::IsMatrix<Eigen::Vector2d>);
    static_assert(math::IsMatrix<Eigen::Vector3d>);
    static_assert(!math::IsMatrix<Eigen::VectorXd&>);
    static_assert(math::IsMatrix<Eigen::MatrixXd>);
    static_assert(math::IsMatrix<Eigen::Matrix<double, 1, 1>>);
    static_assert(math::IsMatrix<Eigen::Matrix2d>);
    static_assert(math::IsMatrix<Eigen::Matrix3d>);
    static_assert(!math::IsMatrix<Eigen::MatrixXd&>);
}

TEST(traits, IsVector) {
    static_assert(!math::IsVector<double>);
    static_assert(!math::IsVector<double&>);
    static_assert(math::IsVector<Eigen::VectorXd>);
    static_assert(math::IsVector<Eigen::Vector<double, 1>>);
    static_assert(math::IsVector<Eigen::Vector2d>);
    static_assert(math::IsVector<Eigen::Vector3d>);
    static_assert(!math::IsVector<Eigen::VectorXd&>);
    static_assert(!math::IsVector<Eigen::MatrixXd>);
    static_assert(math::IsVector<Eigen::Matrix<double, 1, 1>>);
    static_assert(!math::IsVector<Eigen::Matrix2d>);
    static_assert(!math::IsVector<Eigen::Matrix3d>);
    static_assert(!math::IsVector<Eigen::MatrixXd&>);
}

TEST(traits, remove_map) {
    EXPECT_TRUE((std::is_same_v<math::remove_map_t<double>, double>));
    EXPECT_TRUE((std::is_same_v<math::remove_map_t<Eigen::Vector3d>, Eigen::Vector3d>));
    EXPECT_TRUE((std::is_same_v<math::remove_map_t<Eigen::VectorXd>, Eigen::VectorXd>));
    EXPECT_TRUE((std::is_same_v<math::remove_map_t<Eigen::Matrix3d>, Eigen::Matrix3d>));
    EXPECT_TRUE((std::is_same_v<math::remove_map_t<Eigen::MatrixXd>, Eigen::MatrixXd>));
    EXPECT_TRUE((std::is_same_v<math::remove_map_t<Eigen::Map<Eigen::Vector3d>>, Eigen::Vector3d>));
    EXPECT_TRUE((std::is_same_v<math::remove_map_t<Eigen::Map<Eigen::VectorXd>>, Eigen::VectorXd>));
    EXPECT_TRUE((std::is_same_v<math::remove_map_t<Eigen::Map<Eigen::Matrix3d>>, Eigen::Matrix3d>));
    EXPECT_TRUE((std::is_same_v<math::remove_map_t<Eigen::Map<Eigen::MatrixXd>>, Eigen::MatrixXd>));
    EXPECT_TRUE((std::is_same_v<math::remove_map_t<Eigen::Map<const Eigen::Vector3d>>, const Eigen::Vector3d>));
    EXPECT_TRUE((std::is_same_v<math::remove_map_t<Eigen::Map<const Eigen::VectorXd>>, const Eigen::VectorXd>));
    EXPECT_TRUE((std::is_same_v<math::remove_map_t<Eigen::Map<const Eigen::Matrix3d>>, const Eigen::Matrix3d>));
    EXPECT_TRUE((std::is_same_v<math::remove_map_t<Eigen::Map<const Eigen::MatrixXd>>, const Eigen::MatrixXd>));
}

TEST(traits, remove_ref_or_map) {
    EXPECT_TRUE((std::is_same_v<math::remove_ref_or_map_t<double>, double>));
    EXPECT_TRUE((std::is_same_v<math::remove_ref_or_map_t<Eigen::Vector3d>, Eigen::Vector3d>));
    EXPECT_TRUE((std::is_same_v<math::remove_ref_or_map_t<Eigen::VectorXd>, Eigen::VectorXd>));
    EXPECT_TRUE((std::is_same_v<math::remove_ref_or_map_t<Eigen::Matrix3d>, Eigen::Matrix3d>));
    EXPECT_TRUE((std::is_same_v<math::remove_ref_or_map_t<Eigen::MatrixXd>, Eigen::MatrixXd>));
    EXPECT_TRUE((std::is_same_v<math::remove_ref_or_map_t<double&>, double>));
    EXPECT_TRUE((std::is_same_v<math::remove_ref_or_map_t<const double&>, const double>));
    EXPECT_TRUE((std::is_same_v<math::remove_ref_or_map_t<Eigen::Vector3d&>, Eigen::Vector3d>));
    EXPECT_TRUE((std::is_same_v<math::remove_ref_or_map_t<Eigen::VectorXd&>, Eigen::VectorXd>));
    EXPECT_TRUE((std::is_same_v<math::remove_ref_or_map_t<Eigen::Matrix3d&>, Eigen::Matrix3d>));
    EXPECT_TRUE((std::is_same_v<math::remove_ref_or_map_t<Eigen::MatrixXd&>, Eigen::MatrixXd>));
    EXPECT_TRUE((std::is_same_v<math::remove_ref_or_map_t<Eigen::Ref<Eigen::Vector3d>>, Eigen::Vector3d>));
    EXPECT_TRUE((std::is_same_v<math::remove_ref_or_map_t<Eigen::Ref<Eigen::VectorXd>>, Eigen::VectorXd>));
    EXPECT_TRUE((std::is_same_v<math::remove_ref_or_map_t<Eigen::Ref<Eigen::Matrix3d>>, Eigen::Matrix3d>));
    EXPECT_TRUE((std::is_same_v<math::remove_ref_or_map_t<Eigen::Ref<Eigen::MatrixXd>>, Eigen::MatrixXd>));
    EXPECT_TRUE((std::is_same_v<math::remove_ref_or_map_t<Eigen::Ref<const Eigen::Vector3d>>, const Eigen::Vector3d>));
    EXPECT_TRUE((std::is_same_v<math::remove_ref_or_map_t<Eigen::Ref<const Eigen::VectorXd>>, const Eigen::VectorXd>));
    EXPECT_TRUE((std::is_same_v<math::remove_ref_or_map_t<Eigen::Ref<const Eigen::Matrix3d>>, const Eigen::Matrix3d>));
    EXPECT_TRUE((std::is_same_v<math::remove_ref_or_map_t<Eigen::Ref<const Eigen::MatrixXd>>, const Eigen::MatrixXd>));
    EXPECT_TRUE((std::is_same_v<math::remove_ref_or_map_t<Eigen::Map<Eigen::Vector3d>>, Eigen::Vector3d>));
    EXPECT_TRUE((std::is_same_v<math::remove_ref_or_map_t<Eigen::Map<Eigen::VectorXd>>, Eigen::VectorXd>));
    EXPECT_TRUE((std::is_same_v<math::remove_ref_or_map_t<Eigen::Map<Eigen::Matrix3d>>, Eigen::Matrix3d>));
    EXPECT_TRUE((std::is_same_v<math::remove_ref_or_map_t<Eigen::Map<Eigen::MatrixXd>>, Eigen::MatrixXd>));
    EXPECT_TRUE((std::is_same_v<math::remove_ref_or_map_t<Eigen::Map<const Eigen::Vector3d>>, const Eigen::Vector3d>));
    EXPECT_TRUE((std::is_same_v<math::remove_ref_or_map_t<Eigen::Map<const Eigen::VectorXd>>, const Eigen::VectorXd>));
    EXPECT_TRUE((std::is_same_v<math::remove_ref_or_map_t<Eigen::Map<const Eigen::Matrix3d>>, const Eigen::Matrix3d>));
    EXPECT_TRUE((std::is_same_v<math::remove_ref_or_map_t<Eigen::Map<const Eigen::MatrixXd>>, const Eigen::MatrixXd>));
}

TEST(traits, is_ref_or_map) {
    EXPECT_FALSE((math::is_ref_or_map_v<double, double>));
    EXPECT_FALSE((math::is_ref_or_map_v<Eigen::Vector3d, Eigen::Vector3d>));
    EXPECT_FALSE((math::is_ref_or_map_v<Eigen::VectorXd, Eigen::VectorXd>));
    EXPECT_FALSE((math::is_ref_or_map_v<Eigen::Matrix3d, Eigen::Matrix3d>));
    EXPECT_FALSE((math::is_ref_or_map_v<Eigen::MatrixXd, Eigen::MatrixXd>));
    // Normal references to Eigen structures are purposely disallowed
    EXPECT_FALSE((math::is_ref_or_map_v<Eigen::Vector3d&, Eigen::Vector3d>));
    EXPECT_FALSE((math::is_ref_or_map_v<Eigen::VectorXd&, Eigen::VectorXd>));
    EXPECT_FALSE((math::is_ref_or_map_v<Eigen::Matrix3d&, Eigen::Matrix3d>));
    EXPECT_FALSE((math::is_ref_or_map_v<Eigen::MatrixXd&, Eigen::MatrixXd>));
    EXPECT_TRUE((math::is_ref_or_map_v<double&, double>));
    EXPECT_TRUE((math::is_ref_or_map_v<Eigen::Ref<Eigen::Vector3d>, Eigen::Vector3d>));
    EXPECT_TRUE((math::is_ref_or_map_v<Eigen::Ref<Eigen::VectorXd>, Eigen::VectorXd>));
    EXPECT_TRUE((math::is_ref_or_map_v<Eigen::Ref<Eigen::Matrix3d>, Eigen::Matrix3d>));
    EXPECT_TRUE((math::is_ref_or_map_v<Eigen::Ref<Eigen::MatrixXd>, Eigen::MatrixXd>));
    EXPECT_TRUE((math::is_ref_or_map_v<Eigen::Map<Eigen::Vector3d>, Eigen::Vector3d>));
    EXPECT_TRUE((math::is_ref_or_map_v<Eigen::Map<Eigen::VectorXd>, Eigen::VectorXd>));
    EXPECT_TRUE((math::is_ref_or_map_v<Eigen::Map<Eigen::Matrix3d>, Eigen::Matrix3d>));
    EXPECT_TRUE((math::is_ref_or_map_v<Eigen::Map<Eigen::MatrixXd>, Eigen::MatrixXd>));
    EXPECT_TRUE((math::is_ref_or_map_v<const double&, const double>));
    EXPECT_TRUE((math::is_ref_or_map_v<Eigen::Ref<const Eigen::Vector3d>, const Eigen::Vector3d>));
    EXPECT_TRUE((math::is_ref_or_map_v<Eigen::Ref<const Eigen::VectorXd>, const Eigen::VectorXd>));
    EXPECT_TRUE((math::is_ref_or_map_v<Eigen::Ref<const Eigen::Matrix3d>, const Eigen::Matrix3d>));
    EXPECT_TRUE((math::is_ref_or_map_v<Eigen::Ref<const Eigen::MatrixXd>, const Eigen::MatrixXd>));
    EXPECT_TRUE((math::is_ref_or_map_v<Eigen::Map<const Eigen::Vector3d>, const Eigen::Vector3d>));
    EXPECT_TRUE((math::is_ref_or_map_v<Eigen::Map<const Eigen::VectorXd>, const Eigen::VectorXd>));
    EXPECT_TRUE((math::is_ref_or_map_v<Eigen::Map<const Eigen::Matrix3d>, const Eigen::Matrix3d>));
    EXPECT_TRUE((math::is_ref_or_map_v<Eigen::Map<const Eigen::MatrixXd>, const Eigen::MatrixXd>));
    // Check compilation of concept
    static_assert(math::IsRefOrMap<double&, double>);
    static_assert(math::IsRefOrMap<Eigen::Ref<Eigen::Vector3d>, Eigen::Vector3d>);
    static_assert(math::IsRefOrMap<Eigen::Ref<Eigen::VectorXd>, Eigen::VectorXd>);
    static_assert(math::IsRefOrMap<Eigen::Ref<Eigen::Matrix3d>, Eigen::Matrix3d>);
    static_assert(math::IsRefOrMap<Eigen::Ref<Eigen::MatrixXd>, Eigen::MatrixXd>);
    static_assert(math::IsRefOrMap<Eigen::Map<Eigen::Vector3d>, Eigen::Vector3d>);
    static_assert(math::IsRefOrMap<Eigen::Map<Eigen::VectorXd>, Eigen::VectorXd>);
    static_assert(math::IsRefOrMap<Eigen::Map<Eigen::Matrix3d>, Eigen::Matrix3d>);
    static_assert(math::IsRefOrMap<Eigen::Map<Eigen::MatrixXd>, Eigen::MatrixXd>);
}

TEST(traits, is_same_or_const_map) {
    EXPECT_TRUE((math::is_same_or_const_map_v<double, double>));
    EXPECT_TRUE((math::is_same_or_const_map_v<Eigen::Vector3d, Eigen::Vector3d>));
    EXPECT_TRUE((math::is_same_or_const_map_v<Eigen::VectorXd, Eigen::VectorXd>));
    EXPECT_TRUE((math::is_same_or_const_map_v<Eigen::Matrix3d, Eigen::Matrix3d>));
    EXPECT_TRUE((math::is_same_or_const_map_v<Eigen::MatrixXd, Eigen::MatrixXd>));
    EXPECT_TRUE((math::is_same_or_const_map_v<Eigen::Map<const Eigen::Vector3d>, Eigen::Vector3d>));
    EXPECT_TRUE((math::is_same_or_const_map_v<Eigen::Map<const Eigen::VectorXd>, Eigen::VectorXd>));
    EXPECT_TRUE((math::is_same_or_const_map_v<Eigen::Map<const Eigen::Matrix3d>, Eigen::Matrix3d>));
    EXPECT_TRUE((math::is_same_or_const_map_v<Eigen::Map<const Eigen::MatrixXd>, Eigen::MatrixXd>));
    // Check compilation of concept
    static_assert(math::IsSameOrConstMap<double, double>);
    static_assert(math::IsSameOrConstMap<Eigen::Vector3d, Eigen::Vector3d>);
    static_assert(math::IsSameOrConstMap<Eigen::VectorXd, Eigen::VectorXd>);
    static_assert(math::IsSameOrConstMap<Eigen::Matrix3d, Eigen::Matrix3d>);
    static_assert(math::IsSameOrConstMap<Eigen::MatrixXd, Eigen::MatrixXd>);
    static_assert(math::IsSameOrConstMap<Eigen::Map<const Eigen::Vector3d>, Eigen::Vector3d>);
    static_assert(math::IsSameOrConstMap<Eigen::Map<const Eigen::VectorXd>, Eigen::VectorXd>);
    static_assert(math::IsSameOrConstMap<Eigen::Map<const Eigen::Matrix3d>, Eigen::Matrix3d>);
    static_assert(math::IsSameOrConstMap<Eigen::Map<const Eigen::MatrixXd>, Eigen::MatrixXd>);
}

TEST(traits, is_same_or_any_map) {
    EXPECT_TRUE((math::is_same_or_any_map_v<double, double>));
    EXPECT_TRUE((math::is_same_or_any_map_v<Eigen::Vector3d, Eigen::Vector3d>));
    EXPECT_TRUE((math::is_same_or_any_map_v<Eigen::VectorXd, Eigen::VectorXd>));
    EXPECT_TRUE((math::is_same_or_any_map_v<Eigen::Matrix3d, Eigen::Matrix3d>));
    EXPECT_TRUE((math::is_same_or_any_map_v<Eigen::MatrixXd, Eigen::MatrixXd>));
    EXPECT_TRUE((math::is_same_or_any_map_v<Eigen::Map<Eigen::Vector3d>, Eigen::Vector3d>));
    EXPECT_TRUE((math::is_same_or_any_map_v<Eigen::Map<Eigen::VectorXd>, Eigen::VectorXd>));
    EXPECT_TRUE((math::is_same_or_any_map_v<Eigen::Map<Eigen::Matrix3d>, Eigen::Matrix3d>));
    EXPECT_TRUE((math::is_same_or_any_map_v<Eigen::Map<Eigen::MatrixXd>, Eigen::MatrixXd>));
    EXPECT_TRUE((math::is_same_or_any_map_v<Eigen::Map<const Eigen::Vector3d>, Eigen::Vector3d>));
    EXPECT_TRUE((math::is_same_or_any_map_v<Eigen::Map<const Eigen::VectorXd>, Eigen::VectorXd>));
    EXPECT_TRUE((math::is_same_or_any_map_v<Eigen::Map<const Eigen::Matrix3d>, Eigen::Matrix3d>));
    EXPECT_TRUE((math::is_same_or_any_map_v<Eigen::Map<const Eigen::MatrixXd>, Eigen::MatrixXd>));
    // Check compilation of concept
    static_assert(math::IsSameOrAnyMap<double, double>);
    static_assert(math::IsSameOrAnyMap<Eigen::Vector3d, Eigen::Vector3d>);
    static_assert(math::IsSameOrAnyMap<Eigen::VectorXd, Eigen::VectorXd>);
    static_assert(math::IsSameOrAnyMap<Eigen::Matrix3d, Eigen::Matrix3d>);
    static_assert(math::IsSameOrAnyMap<Eigen::MatrixXd, Eigen::MatrixXd>);
    static_assert(math::IsSameOrAnyMap<Eigen::Map<Eigen::Vector3d>, Eigen::Vector3d>);
    static_assert(math::IsSameOrAnyMap<Eigen::Map<Eigen::VectorXd>, Eigen::VectorXd>);
    static_assert(math::IsSameOrAnyMap<Eigen::Map<Eigen::Matrix3d>, Eigen::Matrix3d>);
    static_assert(math::IsSameOrAnyMap<Eigen::Map<Eigen::MatrixXd>, Eigen::MatrixXd>);
    static_assert(math::IsSameOrAnyMap<Eigen::Map<const Eigen::Vector3d>, Eigen::Vector3d>);
    static_assert(math::IsSameOrAnyMap<Eigen::Map<const Eigen::VectorXd>, Eigen::VectorXd>);
    static_assert(math::IsSameOrAnyMap<Eigen::Map<const Eigen::Matrix3d>, Eigen::Matrix3d>);
    static_assert(math::IsSameOrAnyMap<Eigen::Map<const Eigen::MatrixXd>, Eigen::MatrixXd>);
}
