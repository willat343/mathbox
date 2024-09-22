#ifndef MATHBOX_IMPL_NORMAL_HPP
#define MATHBOX_IMPL_NORMAL_HPP

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include "mathbox/normal.hpp"

namespace math {

template<typename Scalar, int Size>
GrvGenerator<Scalar, Size>::GrvGenerator(const Vector& mean, const Matrix& covariance,
        const LLTDecompositionMethod decomposition_method, const unsigned int seed)
    : generator(seed) {
    set_mean_covariance(mean, covariance, decomposition_method);
}

template<typename Scalar, int Size>
auto GrvGenerator<Scalar, Size>::compute_covariance() const -> Matrix {
    return transform() * transform().transpose();
}

template<typename Scalar, int Size>
auto GrvGenerator<Scalar, Size>::mean() const -> const Vector& {
    return mean_;
}

template<typename Scalar, int Size>
auto GrvGenerator<Scalar, Size>::operator()() -> Vector {
    // Note that passing the size to the vector is redundant in the fixed-size case.
    return mean() +
           transform() * Vector(mean_.size())
                                 .unaryExpr([&generator = generator, &normal_distribution = normal_distribution](
                                                    auto) { return normal_distribution(generator); });
}

template<typename Scalar, int Size>
inline void GrvGenerator<Scalar, Size>::set_covariance(const Matrix& covariance,
        const LLTDecompositionMethod decomposition_method) {
    // Check sizes in the dynamic case.
    if constexpr (Size == Eigen::Dynamic) {
        const int mean_size = mean_.size();
        const int covariance_rows = covariance.rows();
        const int covariance_cols = covariance.cols();
        if (covariance_rows != covariance_cols) {
            throw std::runtime_error("Covariance matrix must be square but was " + std::to_string(covariance_rows) +
                                     " X " + std::to_string(covariance_cols) + ".");
        } else if (mean_size != covariance_rows) {
            throw std::runtime_error("Mean vector (" + std::to_string(mean_size) + " X 1) and covariance matrix (" +
                                     std::to_string(covariance_rows) + " X " + std::to_string(covariance_cols) +
                                     ") must have the same dimensionality.");
        }
    }
    set_covariance_impl(covariance, decomposition_method);
}

template<typename Scalar, int Size>
inline void GrvGenerator<Scalar, Size>::set_mean(const Vector& mean) {
    // Check sizes in the dynamic case.
    if constexpr (Size == Eigen::Dynamic) {
        const int mean_size = mean.size();
        const int transform_rows = transform_.rows();
        if (mean_size != transform_rows) {
            throw std::runtime_error("Mean vector (" + std::to_string(mean_size) +
                                     " X 1) must have the same dimensionality as the covariance matrix (" +
                                     std::to_string(transform_rows) + " X " + std::to_string(transform_rows) + ").");
        }
    }
    set_mean_impl(mean);
}

template<typename Scalar, int Size>
inline void GrvGenerator<Scalar, Size>::set_mean_covariance(const Vector& mean, const Matrix& covariance,
        const LLTDecompositionMethod decomposition_method) {
    // Check sizes in the dynamic case.
    if constexpr (Size == Eigen::Dynamic) {
        const int mean_size = mean.size();
        const int covariance_rows = covariance.rows();
        const int covariance_cols = covariance.cols();
        if (covariance_rows != covariance_cols) {
            throw std::runtime_error("Covariance matrix must be square but was " + std::to_string(covariance_rows) +
                                     " X " + std::to_string(covariance_cols) + ".");
        } else if (mean_size != covariance_rows) {
            throw std::runtime_error("Mean vector (" + std::to_string(mean_size) + " X 1) and covariance matrix (" +
                                     std::to_string(covariance_rows) + " X " + std::to_string(covariance_cols) +
                                     ") must have the same dimensionality.");
        }
    }
    set_mean_impl(mean);
    set_covariance_impl(covariance, decomposition_method);
}

template<typename Scalar, int Size>
void GrvGenerator<Scalar, Size>::set_seed(const unsigned int seed) {
    generator = std::mt19937(seed);
}

template<typename Scalar, int Size>
auto GrvGenerator<Scalar, Size>::transform() const -> const Matrix& {
    return transform_;
}

template<typename Scalar, int Size>
inline void GrvGenerator<Scalar, Size>::set_mean_impl(const Vector& mean) {
    mean_ = mean;
}

template<typename Scalar, int Size>
inline void GrvGenerator<Scalar, Size>::set_covariance_impl(const Matrix& covariance,
        const LLTDecompositionMethod decomposition_method) {
    transform_ = LLT(covariance, decomposition_method);
}

}

#endif
