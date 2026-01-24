#ifndef MATHBOX_MATRIX_DIAGNOSTICS_HPP
#define MATHBOX_MATRIX_DIAGNOSTICS_HPP

#include <Eigen/Core>
#include <optional>
#include <string>
#include <type_traits>

#include "mathbox/traits.hpp"

namespace math {

/**
 * @brief Check the computation info and throw an exception if computation was not successful.
 *
 * @param info
 */
void check_computation_info(const Eigen::ComputationInfo info);

/**
 * @brief Produce diagnostic string about the type and size of a matrix. The type is one of:
 *  - Z: integer (additional "^+" symbol indicates the type is unsigned)
 *  - C: complex
 *  - R: real
 *  - 0 (special case): zero matrix
 *  - I (special case): identity matrix
 *
 * The size is appended in form "_NxM".
 *
 * @tparam Derived
 * @param m
 * @return std::string
 */
template<typename Derived>
std::string matrix_type_and_size(const Eigen::MatrixBase<Derived>& m);

/**
 * @brief Produce diagnostic information about the structure of a matrix.
 *
 * @tparam Derived
 * @param m
 * @param row_block_sizes
 * @param column_block_sizes
 * @param row_block_names
 * @param column_block_names
 * @param print_block_norm
 * @return std::string
 */
template<IsMatrix Derived>
std::string structure_diagnostics(const Eigen::MatrixBase<Derived>& m, const std::vector<std::size_t>& row_block_sizes,
        const std::vector<std::size_t>& column_block_sizes,
        const std::optional<std::vector<std::string>>& row_block_names = std::nullopt,
        const std::optional<std::vector<std::string>>& column_block_names = std::nullopt,
        const std::optional<std::string>& matrix_name = std::nullopt, const bool print_block_norm = false);

/**
 * @brief Overload for vector types
 *
 * @tparam Derived
 * @param v
 * @param row_block_sizes
 * @param row_block_names
 * @param vector_name
 * @param print_block_norm
 * @return std::string
 */
template<IsVector Derived>
std::string structure_diagnostics(const Eigen::MatrixBase<Derived>& v, const std::vector<std::size_t>& row_block_sizes,
        const std::optional<std::vector<std::string>>& row_block_names = std::nullopt,
        const std::optional<std::string>& vector_name = std::nullopt, const bool print_block_norm = false);

template<typename Scalar>
struct SpectralStructureSummary {
    template<IsVector Derived>
        requires(std::is_same_v<typename Derived::Scalar, Scalar>)
    SpectralStructureSummary(const Eigen::MatrixBase<Derived>& eigenvalues);

    Scalar min_eigenvalue;
    Scalar max_eigenvalue;
    Scalar nullspace_threshold;
    Scalar weak_threshold;
    std::size_t num_nullspace_eigenvalues;
    std::size_t num_weak_eigenvalues;
    std::size_t num_strong_eigenvalues;
};

/**
 * @brief Produce diagnostic information about the spectral structure of a square matrix given its eigen decomposition.
 *
 * Eigenvector directions are categorised as nullspace, weak, or strong, and their connection with specific elements (or
 * optionally blocks) is reported.
 *
 * @tparam EigenvaluesDerived vector type
 * @tparam EigenvectorsDerived matrix type
 * @param eigenvalues eigenvalues
 * @param eigenvectors eigenvectors, which must be square equal in dimensions to the eigenvalues
 * @param contribution_threshold threshold for an eigenvector element or block norm to be considered a contributor which
 * must be in range [0.0, 1.0], e.g. values in [0.01, 0.1] are typical
 * @param block_sizes optional semantic block sizes, which must sum to the total size (i.e., number of eigenvectors) and
 * cannot have zero elements
 * @param block_names optional semantic block names, which if defined requires that `block_sizes` be defined and must be
 * equal in size to block sizes
 * @param matrix_name optional matrix name
 * @param exclude_nullspaces if true, exclude nullspaces from reported diagnostics
 * @param exclude_weak if true, exclude weak directions from reported diagnostics
 * @param exclude_strong if true, exclude strong directions from reported diagnostics
 * @return std::string diagnostic information about the nullspace structure
 */
template<IsVector EigenvaluesDerived, IsMatrix EigenvectorsDerived>
    requires(std::is_same_v<typename EigenvaluesDerived::Scalar, typename EigenvectorsDerived::Scalar>)
std::string spectral_structure_diagnostics(const Eigen::MatrixBase<EigenvaluesDerived>& eigenvalues,
        const Eigen::MatrixBase<EigenvectorsDerived>& eigenvectors,
        const typename EigenvaluesDerived::Scalar contribution_threshold,
        const std::optional<std::vector<std::size_t>>& block_sizes = std::nullopt,
        const std::optional<std::vector<std::string>>& block_names = std::nullopt,
        const std::optional<std::string>& matrix_name = std::nullopt, const bool exclude_nullspaces = false,
        const bool exclude_weak = false, const bool exclude_strong = false);

}

#include "mathbox/impl/matrix_diagnostics.hpp"

#endif
