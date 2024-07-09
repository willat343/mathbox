#ifndef MATHBOX_STIFFNESS_HPP
#define MATHBOX_STIFFNESS_HPP

#include <Eigen/Core>

#include "mathbox/traits.hpp"

namespace math {

/**
 * @brief Compute the stiffness matrix \f$\mathbf{A}\f$ from the covariance matrix \f$\boldsymbol{\Sigma}\f$ such that
 * \f$\boldsymbol{\Lambda} = \boldsymbol{\Sigma}^{-1} = \mathbf{A}^T\mathbf{A}\f$, with order preserved so that
 * computing the cost of a residual with this stiffness matrix as \f$\mathbf{A} \mathbf{e}\f$ for Mahalanobis norm
 * \f$\Vert \mathbf{e} \Vert_{\boldsymbol{\Sigma}}^2 = \mathbf{e}^T \boldsymbol{\Sigma}^{-1} \mathbf{e} =
 * \mathbf{e}^T \mathbf{A}^T \mathbf{A} \mathbf{e} = (\mathbf{A} \mathbf{e})^T (\mathbf{A} \mathbf{e})\f$
 * preserves the order of the terms in \f$\mathbf{e}\f$.
 *
 * @tparam Derived
 * @param covariance
 * @return Derived
 */
template<typename Derived>
Derived stiffness_from_covariance(const Eigen::MatrixBase<Derived>& covariance);

/**
 * @brief Alternative function which allows calls without creating temporaries, however requires explicit template
 * parameter specification.
 *
 * @tparam Rows
 * @tparam Scalar
 * @param covariance
 * @return Eigen::Matrix<Scalar, Rows, Rows>
 */
template<int Rows, typename Scalar = double>
Eigen::Matrix<Scalar, Rows, Rows> stiffness_from_covariance(
        const Eigen::Ref<const Eigen::Matrix<Scalar, Rows, Rows>>& covariance);

/**
 * @brief Compute the stiffness matrix \f$\mathbf{A}\f$ from the information matrix \f$\boldsymbol{\Lambda}\f$ such that
 * \f$\boldsymbol{\Lambda} = \boldsymbol{\Sigma}^{-1} = \mathbf{A}^T\mathbf{A}\f$, with order preserved so that
 * computing the cost of a residual with this stiffness matrix as \f$\mathbf{A} \mathbf{e}\f$ for Mahalanobis norm
 * \f$\Vert \mathbf{e} \Vert_{\boldsymbol{\Sigma}}^2 = \mathbf{e}^T \boldsymbol{\Sigma}^{-1} \mathbf{e} =
 * \mathbf{e}^T \mathbf{A}^T \mathbf{A} \mathbf{e} = (\mathbf{A} \mathbf{e})^T (\mathbf{A} \mathbf{e})\f$
 * preserves the order of the terms in \f$\mathbf{e}\f$.
 *
 * @tparam Derived
 * @param information
 * @return Derived
 */
template<typename Derived>
Derived stiffness_from_information(const Eigen::MatrixBase<Derived>& information);

/**
 * @brief Alternative function which allows calls without creating temporaries, however requires explicit template
 * parameter specification.
 *
 * @tparam Rows
 * @tparam Scalar
 * @param information
 * @return Eigen::Matrix<Scalar, Rows, Rows>
 */
template<int Rows, typename Scalar = double>
Eigen::Matrix<Scalar, Rows, Rows> stiffness_from_information(
        const Eigen::Ref<const Eigen::Matrix<Scalar, Rows, Rows>>& information);

/**
 * @brief Compute the stiffness matrix from a vector of standard deviations (sigmas), the result being a diagonal matrix
 * with the inverse of these sigmas along the diagonal.
 *
 * @tparam Derived
 * @param sigmas
 * @return Eigen::Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, Derived::RowsAtCompileTime>
 */
template<typename Derived>
Eigen::Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, Derived::RowsAtCompileTime> stiffness_from_sigmas(
        const Eigen::MatrixBase<Derived>& sigmas);

/**
 * @brief Alternative function which allows calls without creating temporaries, however requires explicit template
 * parameter specification.
 *
 * @tparam Rows
 * @tparam Scalar
 * @param sigmas
 * @return Eigen::Matrix<Scalar, Rows, Rows>
 */
template<int Rows, typename Scalar = double>
Eigen::Matrix<Scalar, Rows, Rows> stiffness_from_sigmas(const Eigen::Ref<const Eigen::Matrix<Scalar, Rows, 1>>& sigmas);

/**
 * @brief Like stiffness_from_sigmas, except where all the standard deviations are the same. The size must therefore be
 * provided as a template parameter.
 *
 * @tparam Rows
 * @tparam Scalar
 * @param sigma
 * @return Eigen::Matrix<Scalar, Rows, Rows>
 */
template<int Rows, typename Scalar = double>
Eigen::Matrix<Scalar, Rows, Rows> stiffness_from_sigma(const Scalar sigma);

/**
 * @brief Compute stiffness from a single sigma value.
 *
 * @tparam Scalar
 * @param sigma
 * @return Scalar
 */
template<typename Scalar = double>
Scalar stiffness_from_sigma(const Scalar sigma);

/**
 * @brief Compute stiffness from variances.
 *
 * @tparam Derived
 * @param variances
 * @return Eigen::Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, Derived::RowsAtCompileTime>
 */
template<typename Derived>
Eigen::Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, Derived::RowsAtCompileTime>
stiffness_from_variances(const Eigen::MatrixBase<Derived>& variances);

/**
 * @brief Alternative function which allows calls without creating temporaries, however requires explicit template
 * parameter specification.
 *
 * @tparam Rows
 * @tparam Scalar
 * @param variances
 * @return Eigen::Matrix<Scalar, Rows, Rows>
 */
template<int Rows, typename Scalar = double>
Eigen::Matrix<Scalar, Rows, Rows> stiffness_from_variances(
        const Eigen::Ref<const Eigen::Matrix<Scalar, Rows, 1>>& variances);

/**
 * @brief Like stiffness_from_variances, except where all the variances are the same. The size must therefore be
 * provided as a template parameter.
 *
 * @tparam Rows
 * @tparam Scalar
 * @param variance
 * @return Eigen::Matrix<Scalar, Rows, Rows>
 */
template<int Rows, typename Scalar = double>
Eigen::Matrix<Scalar, Rows, Rows> stiffness_from_variance(const Scalar variance);

/**
 * @brief Compute stiffness from a single variance value.
 *
 * @tparam Scalar
 * @param variance
 * @return Scalar
 */
template<typename Scalar = double>
Scalar stiffness_from_variance(const Scalar variance);

}

#include "mathbox/impl/stiffness.hpp"

#endif
