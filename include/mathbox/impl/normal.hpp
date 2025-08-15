#ifndef MATHBOX_IMPL_NORMAL_HPP
#define MATHBOX_IMPL_NORMAL_HPP

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <cppbox/exceptions.hpp>

#include "mathbox/normal.hpp"

namespace math {

template<typename Scalar_, int Size_>
GrvGenerator<Scalar_, Size_>::GrvGenerator(const Vector& mean, const Matrix& covariance,
        const LLTDecompositionMethod decomposition_method, const unsigned int seed)
    : generator(seed) {
    set_mean_covariance(mean, covariance, decomposition_method);
}

template<typename Scalar_, int Size_>
auto GrvGenerator<Scalar_, Size_>::compute_covariance() const -> Matrix {
    return transform() * transform().transpose();
}

template<typename Scalar_, int Size_>
auto GrvGenerator<Scalar_, Size_>::mean() const -> const Vector& {
    return mean_;
}

template<typename Scalar_, int Size_>
auto GrvGenerator<Scalar_, Size_>::operator()() -> Vector {
    // Note that passing the size to the vector is redundant in the fixed-size case.
    return mean() +
           transform() * Vector(mean_.size())
                                 .unaryExpr([&generator = generator, &normal_distribution = normal_distribution](
                                                    auto) { return normal_distribution(generator); });
}

template<typename Scalar_, int Size_>
inline void GrvGenerator<Scalar_, Size_>::set_covariance(const Matrix& covariance,
        const LLTDecompositionMethod decomposition_method) {
    // Check sizes in the dynamic case.
    if constexpr (Size == Eigen::Dynamic) {
        const int mean_size = mean_.size();
        const int covariance_rows = covariance.rows();
        const int covariance_cols = covariance.cols();
        throw_if(covariance_rows != covariance_cols, "Covariance matrix must be square but was " +
                                                             std::to_string(covariance_rows) + " X " +
                                                             std::to_string(covariance_cols) + ".");
        throw_if(mean_size != covariance_rows,
                "Mean vector (" + std::to_string(mean_size) + " X 1) and covariance matrix (" +
                        std::to_string(covariance_rows) + " X " + std::to_string(covariance_cols) +
                        ") must have the same dimensionality.");
    }
    set_covariance_impl(covariance, decomposition_method);
}

template<typename Scalar_, int Size_>
inline void GrvGenerator<Scalar_, Size_>::set_mean(const Vector& mean) {
    // Check sizes in the dynamic case.
    if constexpr (Size == Eigen::Dynamic) {
        const int mean_size = mean.size();
        const int transform_rows = transform_.rows();
        throw_if(mean_size != transform_rows,
                "Mean vector (" + std::to_string(mean_size) +
                        " X 1) must have the same dimensionality as the covariance matrix (" +
                        std::to_string(transform_rows) + " X " + std::to_string(transform_rows) + ").");
    }
    set_mean_impl(mean);
}

template<typename Scalar_, int Size_>
inline void GrvGenerator<Scalar_, Size_>::set_mean_covariance(const Vector& mean, const Matrix& covariance,
        const LLTDecompositionMethod decomposition_method) {
    // Check sizes in the dynamic case.
    if constexpr (Size == Eigen::Dynamic) {
        const int mean_size = mean.size();
        const int covariance_rows = covariance.rows();
        const int covariance_cols = covariance.cols();
        throw_if(covariance_rows != covariance_cols, "Covariance matrix must be square but was " +
                                                             std::to_string(covariance_rows) + " X " +
                                                             std::to_string(covariance_cols) + ".");
        throw_if(mean_size != covariance_rows,
                "Mean vector (" + std::to_string(mean_size) + " X 1) and covariance matrix (" +
                        std::to_string(covariance_rows) + " X " + std::to_string(covariance_cols) +
                        ") must have the same dimensionality.");
    }
    set_mean_impl(mean);
    set_covariance_impl(covariance, decomposition_method);
}

template<typename Scalar_, int Size_>
void GrvGenerator<Scalar_, Size_>::set_seed(const unsigned int seed) {
    generator = std::mt19937(seed);
}

template<typename Scalar_, int Size_>
auto GrvGenerator<Scalar_, Size_>::transform() const -> const Matrix& {
    return transform_;
}

template<typename Scalar_, int Size_>
inline void GrvGenerator<Scalar_, Size_>::set_mean_impl(const Vector& mean) {
    mean_ = mean;
}

template<typename Scalar_, int Size_>
inline void GrvGenerator<Scalar_, Size_>::set_covariance_impl(const Matrix& covariance,
        const LLTDecompositionMethod decomposition_method) {
    transform_ = LLT(covariance, decomposition_method);
}

}

#endif
