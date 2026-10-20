#ifndef TEST_IMPL_HELPERS_HPP
#define TEST_IMPL_HELPERS_HPP

#include <cmath>

#include "test/helpers.hpp"

inline Eigen::MatrixXd well_conditioned_spd(const int size) {
    Eigen::MatrixXd M(size, size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            M(i, j) = std::sin(static_cast<double>(size * i + j + 1));
        }
    }
    return M * M.transpose() + static_cast<double>(size) * Eigen::MatrixXd::Identity(size, size);
}

inline Eigen::VectorXd deterministic_vector(const int size) {
    Eigen::VectorXd v(size);
    for (int i = 0; i < size; ++i) {
        v[i] = std::cos(static_cast<double>(i + 1));
    }
    return v;
}

inline Eigen::MatrixXd deterministic_matrix(const int rows, const int cols) {
    Eigen::MatrixXd M(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            M(i, j) = std::sin(static_cast<double>((i + 1) * (j + 1)));
        }
    }
    return M;
}

template<math::IsEuclideanType EuclideanType>
EuclideanType zero_instance() {
    if constexpr (math::IsMatrix<EuclideanType>) {
        // Set any dynamic-sized dimensions to 1
        return math::MathTypeTraits<EuclideanType>::zero(
                EuclideanType::RowsAtCompileTime == Eigen::Dynamic ? 1 : EuclideanType::RowsAtCompileTime,
                EuclideanType::ColsAtCompileTime == Eigen::Dynamic ? 1 : EuclideanType::ColsAtCompileTime);
    } else {
        return math::MathTypeTraits<EuclideanType>::zero();
    }
}

#endif
