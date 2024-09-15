#ifndef MATHBOX_NORMAL_HPP
#define MATHBOX_NORMAL_HPP

#include <Eigen/Core>
#include <random>
#include <type_traits>
#include <vector>

#include "mathbox/decompose.hpp"

namespace math {

/**
 * @brief A Gaussian Random Variable (GRV) generator, given its mean and covariance.
 *
 * @tparam Scalar_ floating point type
 * @tparam Size_ vector dimension, which may be Eigen::Dynamic
 */
template<typename Scalar_, int Size_>
class GrvGenerator {
public:
    using Scalar = Scalar_;
    static constexpr int Size = Size_;
    using DiagonalMatrix = Eigen::DiagonalMatrix<Scalar, Size>;
    using Matrix = Eigen::Matrix<Scalar, Size, Size>;
    using Vector = Eigen::Matrix<Scalar, Size, 1>;

    // Type requirements
    static_assert(
            std::is_same_v<Scalar, float> || std::is_same_v<Scalar, double> || std::is_same_v<Scalar, long double>,
            "Scalar must be a one of float, double or long double because any other types are undefined for "
            "std::normal_distribution.");
    static_assert(Size > 0 || Size == Eigen::Dynamic, "Size must be greater than 0 or Eigen::Dynamic.");

    /**
     * @brief Construct a Gaussian Random Variable (GRV) generator.
     *
     * @param mean mean vector
     * @param covariance covariance matrix
     * @param decomposition_method method to decompose the covariance matrix
     * @param seed seed for the PRNG
     */
    explicit GrvGenerator(const Vector& mean, const Matrix& covariance,
            const LLTDecompositionMethod decomposition_method = LLTDecompositionMethod::ROBUST_CHOLESKY,
            const unsigned int seed = std::random_device{}());

    /**
     * @brief Recover the covariance matrix from its stored decomposition, equivalent to `transform() *
     * transform().transpose()`.
     *
     * @return Matrix
     */
    Matrix compute_covariance() const;

    /**
     * @brief Get the mean vector.
     *
     * @return const Vector&
     */
    const Vector& mean() const;

    /**
     * @brief Generate a new GRV sample \f$\mathbf{x}\f$, where \f$\mathbf{y} = \mathcal{N}(\mathbf{0}, \mathbf{I})\f$,
     * as:
     *
     * \f[
     *      \mathbf{x} = \begin{bmatrix} X_1 \\ X_2 \\ \vdots X_N\end{bmatrix} = \boldsymbol{\mu}_\mathbf{x} +
     *          \mathbf{C} \mathbf{y} \sim \mathcal{N}(\boldsymbol{\mu}_\mathbf{x}, \boldsymbol{\Sigma}_\mathbf{x})
     * \f]
     *
     * @return Vector the new GRV sample \f$\mathbf{x}\f$
     */
    Vector operator()();

    /**
     * @brief Update the generator according to a new covariance matrix.
     *
     * @param covariance
     * @param decomposition_method
     */
    void set_covariance(const Matrix& covariance,
            const LLTDecompositionMethod decomposition_method = LLTDecompositionMethod::ROBUST_CHOLESKY);

    /**
     * @brief Update the generator according to a new mean.
     *
     * @param mean
     */
    void set_mean(const Vector& mean);

    /**
     * @brief Set the generator according to a new mean and covariance.
     *
     * @param mean
     * @param covariance
     * @param decomposition_method
     */
    void set_mean_covariance(const Vector& mean, const Matrix& covariance,
            const LLTDecompositionMethod decomposition_method = LLTDecompositionMethod::ROBUST_CHOLESKY);

    /**
     * @brief Set the random seed of the PRNG.
     *
     * @param seed
     */
    void set_seed(const unsigned int seed);

    /**
     * @brief Get the internal matrix \f$\mathbf{C}\f$ computed from the decomposition of the covariance matrix
     * \f$\boldsymbol{\Sigma}\f$ such that \f$\boldsymbol{\Sigma} = \mathbf{C}\mathbf{C}^T\f$.
     *
     * @return const Matrix&
     */
    const Matrix& transform() const;

private:
    /**
     * @brief Pseudo-Random Number Generator (PRNG) for the normal distribution.
     *
     */
    std::mt19937 generator;

    /**
     * @brief Standard normal distribution.
     *
     */
    std::normal_distribution<Scalar> normal_distribution;

    /**
     * @brief Mean vector.
     *
     */
    Vector mean_;

    /**
     * @brief Matrix obtained from the decomposition of the covariance matrix used to transform a GRV of indepedent
     * standard normal variables.
     *
     */
    Matrix transform_;

    /**
     * @brief Implementation of mean update.
     *
     * @param mean
     */
    void set_mean_impl(const Vector& mean);

    /**
     * @brief Implementation of covariance update.
     *
     * @param covariance
     * @param decomposition_method
     */
    void set_covariance_impl(const Matrix& covariance, const LLTDecompositionMethod decomposition_method);
};

}

#include "mathbox/impl/normal.hpp"

#endif
