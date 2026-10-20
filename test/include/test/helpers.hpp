#ifndef TEST_HELPERS_HPP
#define TEST_HELPERS_HPP

#include <Eigen/Core>

#include "mathbox/traits.hpp"

/**
 * @brief Deterministic symmetric positive definite matrix \f$ M M^T + n I \f$ of size \f$ n \f$, with
 * \f$ M(i, j) = \sin(n i + j + 1) \f$. Its entries are non-dyadic, so that different floating-point evaluation paths
 * round differently and can only be compared approximately. Note that \f$ M \f$ alone has rank 2.
 *
 * @param size
 * @return Eigen::MatrixXd
 */
Eigen::MatrixXd well_conditioned_spd(const int size);

/**
 * @brief Deterministic vector with entries \f$ v(i) = \cos(i + 1) \f$.
 *
 * @param size
 * @return Eigen::VectorXd
 */
Eigen::VectorXd deterministic_vector(const int size);

/**
 * @brief Deterministic full-rank, well-conditioned matrix with entries \f$ M(i, j) = \sin((i + 1) (j + 1)) \f$,
 * unlike \f$ \sin(a i + b j + c) \f$, whose angle sum identity makes any such matrix rank 2.
 *
 * @param rows
 * @param cols
 * @return Eigen::MatrixXd
 */
Eigen::MatrixXd deterministic_matrix(const int rows, const int cols);

/**
 * @brief Zero instance of a Euclidean type, with any dynamic-sized dimensions set to 1.
 *
 * @tparam EuclideanType
 * @return EuclideanType
 */
template<math::IsEuclideanType EuclideanType>
EuclideanType zero_instance();

#include "test/impl/helpers.hpp"

#endif
