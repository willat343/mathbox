#include "mathbox/vector_operations.hpp"

#include <gtest/gtest.h>

TEST(vector_operations, lin_spaced_0) {
    const double start{0.0}, end{4.0}, step{1.0};
    const Eigen::VectorXd result = math::lin_spaced_vector(step, start, end);
    const int expected_size{5};
    const double expected_step{1.0};
    Eigen::VectorXd expected(expected_size);
    for (Eigen::Index i = 0; i < expected_size; ++i) {
        expected[i] = start + static_cast<double>(i) * expected_step;
    }
    EXPECT_EQ(result.size(), expected.size());
    EXPECT_TRUE(result.isApprox(expected));
}

TEST(vector_operations, lin_spaced_1) {
    const double start{1.0}, end{5.1}, step{1.0};
    const Eigen::VectorXd result = math::lin_spaced_vector(step, start, end);
    const int expected_size{5};
    const double expected_step{1.025};
    Eigen::VectorXd expected(expected_size);
    for (Eigen::Index i = 0; i < expected_size; ++i) {
        expected[i] = start + static_cast<double>(i) * expected_step;
    }
    EXPECT_EQ(result.size(), expected.size());
    EXPECT_TRUE(result.isApprox(expected));
}

TEST(vector_operations, lin_spaced_2) {
    const double start{1.3}, end{11.2}, step{2.0};
    const Eigen::VectorXd result = math::lin_spaced_vector(step, start, end);
    const int expected_size{5};
    const double expected_step{2.475};
    Eigen::VectorXd expected(expected_size);
    for (Eigen::Index i = 0; i < expected_size; ++i) {
        expected[i] = start + static_cast<double>(i) * expected_step;
    }
    EXPECT_EQ(result.size(), expected.size());
    EXPECT_TRUE(result.isApprox(expected));
}

TEST(vector_operations, lin_spaced_3) {
    const double start{-1.3}, end{-11.2}, step{-2.0};
    const Eigen::VectorXd result = math::lin_spaced_vector(step, start, end);
    const int expected_size{5};
    const double expected_step{-2.475};
    Eigen::VectorXd expected(expected_size);
    for (Eigen::Index i = 0; i < expected_size; ++i) {
        expected[i] = start + static_cast<double>(i) * expected_step;
    }
    EXPECT_EQ(result.size(), expected.size());
    EXPECT_TRUE(result.isApprox(expected));
}

TEST(vector_operations, lin_spaced_4) {
    const double start{-11.2}, end{-1.3}, step{2.0};
    const Eigen::VectorXd result = math::lin_spaced_vector(step, start, end);
    const int expected_size{5};
    const double expected_step{2.475};
    Eigen::VectorXd expected(expected_size);
    for (Eigen::Index i = 0; i < expected_size; ++i) {
        expected[i] = start + static_cast<double>(i) * expected_step;
    }
    EXPECT_EQ(result.size(), expected.size());
    EXPECT_TRUE(result.isApprox(expected));
}

TEST(vector_operations, range_0) {
    const double start{0.0}, end{4.0}, step{1.0};
    const Eigen::VectorXd result = math::range(step, start, end);
    const int expected_size{5};
    const double expected_step = step;
    Eigen::VectorXd expected(expected_size);
    for (Eigen::Index i = 0; i < expected_size; ++i) {
        expected[i] = start + static_cast<double>(i) * expected_step;
    }
    EXPECT_EQ(result.size(), expected.size());
    EXPECT_TRUE(result.isApprox(expected));
}

TEST(vector_operations, range_1) {
    const double start{1.0}, end{5.1}, step{1.0};
    const Eigen::VectorXd result = math::range(step, start, end);
    const int expected_size{5};
    const double expected_step = step;
    Eigen::VectorXd expected(expected_size);
    for (Eigen::Index i = 0; i < expected_size; ++i) {
        expected[i] = start + static_cast<double>(i) * expected_step;
    }
    EXPECT_EQ(result.size(), expected.size());
    EXPECT_TRUE(result.isApprox(expected));
}

TEST(vector_operations, range_2) {
    const double start{1.3}, end{11.2}, step{2.0};
    const Eigen::VectorXd result = math::range(step, start, end);
    const int expected_size{5};
    const double expected_step = step;
    Eigen::VectorXd expected(expected_size);
    for (Eigen::Index i = 0; i < expected_size; ++i) {
        expected[i] = start + static_cast<double>(i) * expected_step;
    }
    EXPECT_EQ(result.size(), expected.size());
    EXPECT_TRUE(result.isApprox(expected));
}

TEST(vector_operations, range_3) {
    const double start{-1.3}, end{-11.2}, step{-2.0};
    const Eigen::VectorXd result = math::range(step, start, end);
    const int expected_size{5};
    const double expected_step = step;
    Eigen::VectorXd expected(expected_size);
    for (Eigen::Index i = 0; i < expected_size; ++i) {
        expected[i] = start + static_cast<double>(i) * expected_step;
    }
    EXPECT_EQ(result.size(), expected.size());
    EXPECT_TRUE(result.isApprox(expected));
}

TEST(vector_operations, range_4) {
    const double start{-11.2}, end{-1.3}, step{2.0};
    const Eigen::VectorXd result = math::range(step, start, end);
    const int expected_size{5};
    const double expected_step = step;
    Eigen::VectorXd expected(expected_size);
    for (Eigen::Index i = 0; i < expected_size; ++i) {
        expected[i] = start + static_cast<double>(i) * expected_step;
    }
    EXPECT_EQ(result.size(), expected.size());
    EXPECT_TRUE(result.isApprox(expected));
}
