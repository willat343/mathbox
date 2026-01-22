#ifndef MATHBOX_IMPL_MATRIX_DIAGNOSTICS_HPP
#define MATHBOX_IMPL_MATRIX_DIAGNOSTICS_HPP

#include <algorithm>
#include <cppbox/parse.hpp>
#include <numeric>
#include <sstream>
#include <utility>

#include "mathbox/matrix_diagnostics.hpp"

namespace math {

inline void check_computation_info(const Eigen::ComputationInfo info) {
    if (info != Eigen::ComputationInfo::Success) {
        throw_if(info == Eigen::ComputationInfo::NumericalIssue, "NumericalIssue encountered during computation.");
        throw_if(info == Eigen::ComputationInfo::NoConvergence, "NoConvergence encountered during computation.");
        throw_if(info == Eigen::InvalidInput, "InvalidInput encountered during computation.");
        throw_here("Unknown error encountered during computation.");
    }
}

template<typename Derived>
std::string matrix_type_and_size(const Eigen::MatrixBase<Derived>& m) {
    std::stringstream ss;
    if (m.isZero()) {
        ss << "0";
    } else if (m.isIdentity()) {
        ss << "I";
    } else {
        if (Eigen::NumTraits<typename Derived::Scalar>::IsInteger) {
            ss << "Z";
        } else if (Eigen::NumTraits<typename Derived::Scalar>::IsComplex) {
            ss << "C";
        } else {
            ss << "R";
        }
        if (!Eigen::NumTraits<typename Derived::Scalar>::IsSigned) {
            ss << "^+";
        }
    }
    ss << "_" << m.rows() << "x" << m.cols();
    return ss.str();
}

template<IsMatrix Derived>
std::string structure_diagnostics(const Eigen::MatrixBase<Derived>& m, const std::vector<std::size_t>& row_block_sizes,
        const std::vector<std::size_t>& column_block_sizes,
        const std::optional<std::vector<std::string>>& row_block_names,
        const std::optional<std::vector<std::string>>& column_block_names,
        const std::optional<std::string>& matrix_name, const bool print_block_norm) {
    // Check validity of inputs
    assert(!row_block_names || row_block_sizes.size() == row_block_names->size());
    assert(!column_block_names || column_block_sizes.size() == column_block_names->size());
    assert(std::count(row_block_sizes.begin(), row_block_sizes.end(), 0) > 0 ||
            std::accumulate(row_block_sizes.begin(), row_block_sizes.end(), 0) == m.rows());
    assert(std::count(column_block_sizes.begin(), column_block_sizes.end(), 0) > 0 ||
            std::accumulate(column_block_sizes.begin(), column_block_sizes.end(), 0) == m.cols());

    // Generate diagnostic information
    std::stringstream ss;
    if (m.size() > 0) {
        // Precompute min/max string lengths
        const std::size_t max_row_index_width =
                std::to_string(row_block_sizes.empty() ? 0 : row_block_sizes.size() - 1).size();
        const std::size_t min_column_width = 5 + 2 * max_row_index_width + (print_block_norm ? 12 : 0);

        // Row identifiers
        std::vector<std::string> row_identifiers(row_block_sizes.size());
        for (std::size_t i = 0; i < row_block_sizes.size(); ++i) {
            std::stringstream internal_ss;
            internal_ss << std::setw(min_column_width)
                        << (row_block_names && !(*row_block_names)[i].empty() ? (*row_block_names)[i] : "Block") << "("
                        << row_block_sizes[i] << ")";
            row_identifiers[i] = internal_ss.str();
        }
        const std::size_t max_row_identifier_width = cppbox::max_size(row_identifiers);

        // Column identifiers
        std::vector<std::string> column_identifiers(column_block_sizes.size());
        for (std::size_t i = 0; i < column_block_sizes.size(); ++i) {
            std::stringstream internal_ss;
            internal_ss << std::setw(min_column_width)
                        << (column_block_names && !(*column_block_names)[i].empty() ? (*column_block_names)[i]
                                                                                    : "Block") +
                                   "(" + std::to_string(column_block_sizes[i]) + ")";
            column_identifiers[i] = internal_ss.str();
        }

        // Column header
        std::stringstream column_header_ss;
        column_header_ss << std::string(max_row_index_width + 3 + max_row_identifier_width, ' ');
        for (std::size_t i = 0; i < column_block_sizes.size(); ++i) {
            column_header_ss << " | " << column_identifiers[i];
        }
        column_header_ss << " |\n";
        const std::string column_header = column_header_ss.str();

        // Header
        const std::string header_start = "=== " + matrix_name.value_or("Matrix") + " Structure ";
        ss << header_start << std::string(column_header.size() - header_start.size() - 1, '=') << "\n";
        ss << column_header;

        // Information for each block
        for (std::size_t i = 0, r = 0; i < row_block_sizes.size(); ++i, r += row_block_sizes[i]) {
            ss << std::setw(max_row_index_width) << i << " | " << std::setw(max_row_identifier_width)
               << row_identifiers[i];
            for (std::size_t j = 0, c = 0; j < column_block_sizes.size(); ++j, c += column_block_sizes[j]) {
                const Eigen::Ref<const Eigen::MatrixXd>& block =
                        m.block(r, c, row_block_sizes[i], column_block_sizes[j]);
                std::stringstream block_info_ss;
                block_info_ss << matrix_type_and_size(block);
                if (print_block_norm) {
                    block_info_ss << " (" << std::scientific << std::setprecision(3) << block.norm() << ")";
                }
                ss << " | " << std::setw(std::max(min_column_width, column_identifiers[j].size()))
                   << block_info_ss.str();
            }
            ss << " |\n";
        }
        ss << std::string(column_header.size() - 1, '=') << "\n";
    }
    return ss.str();
}

template<IsVector Derived>
inline std::string structure_diagnostics(const Eigen::MatrixBase<Derived>& v,
        const std::vector<std::size_t>& row_block_sizes, const std::optional<std::vector<std::string>>& row_block_names,
        const std::optional<std::string>& vector_name, const bool print_block_norm) {
    return structure_diagnostics(v, row_block_sizes, std::vector<std::size_t>(1, 1), row_block_names, std::nullopt,
            vector_name, print_block_norm);
}

template<IsVector EigenvaluesDerived, IsMatrix EigenvectorsDerived>
    requires(std::is_same_v<typename EigenvaluesDerived::Scalar, typename EigenvectorsDerived::Scalar>)
std::string spectral_structure_diagnostics(const Eigen::MatrixBase<EigenvaluesDerived>& eigenvalues,
        const Eigen::MatrixBase<EigenvectorsDerived>& eigenvectors,
        const typename EigenvaluesDerived::Scalar nullspace_threshold,
        const typename EigenvaluesDerived::Scalar weak_threshold,
        const typename EigenvaluesDerived::Scalar contribution_threshold,
        const std::optional<std::vector<std::size_t>>& block_sizes,
        const std::optional<std::vector<std::string>>& block_names, const std::optional<std::string>& matrix_name,
        const bool exclude_nullspaces, const bool exclude_weak, const bool exclude_strong) {
    using Scalar = typename EigenvaluesDerived::Scalar;
    // Check validity of inputs
    assert(eigenvalues.size() == eigenvectors.rows());
    assert(eigenvectors.rows() == eigenvectors.cols());
    assert(!block_sizes || std::count(block_sizes->begin(), block_sizes->end(), 0) > 0 ||
            std::accumulate(block_sizes->begin(), block_sizes->end(), 0) == eigenvalues.size());
    assert(!block_names || (block_sizes && block_sizes->size() == block_names->size()));
    assert(nullspace_threshold >= static_cast<Scalar>(0));
    assert(weak_threshold >= nullspace_threshold);
    assert(contribution_threshold >= static_cast<Scalar>(0) && contribution_threshold <= static_cast<Scalar>(1));

    // Count eigenvalues in each category
    const int num_nullspaces = (eigenvalues.array() < nullspace_threshold).count();
    const int num_weak = (eigenvalues.array() < weak_threshold).count() - num_nullspaces;
    const int num_strong = eigenvalues.size() - num_weak - num_nullspaces;

    // Generate diagnostic information
    std::stringstream ss;
    if (eigenvalues.size() > 0) {
        // Header
        ss << std::scientific << std::setprecision(3) << "=== " << matrix_name.value_or("Matrix")
           << " Spectral Structure [Max Eigenvalue: " << eigenvalues.maxCoeff()
           << ", Min Eigenvalue: " << eigenvalues.minCoeff() << ", " << num_strong << " Strong"
           << ", " << num_weak << " Weak (Threshold: " << weak_threshold << ")"
           << ", " << num_nullspaces << " Nullspaces (Threshold: " << nullspace_threshold << ")] ===\n";
        const std::size_t header_width = ss.str().size() - 1;
        const std::size_t max_index_width = std::to_string(eigenvalues.size() - 1).size();

        // Information for each eigenvalue
        for (int i = 0; i < eigenvalues.size(); ++i) {
            // Eigenvalue and Eigenvector
            const Scalar eigenvalue = eigenvalues[i];
            const Eigen::Ref<const Eigen::Vector<Scalar, EigenvectorsDerived::RowsAtCompileTime>>& eigenvector =
                    eigenvectors.col(i);

            // Get significant contributors of eigenvector
            std::vector<std::pair<std::string, Scalar>> contributions;
            std::size_t block_index{0};
            for (int j = 0; j < eigenvector.size(); j += block_sizes ? (*block_sizes)[block_index] : 1, ++block_index) {
                const Eigen::Ref<const Eigen::Vector<Scalar, Eigen::Dynamic>> eigenvector_block =
                        eigenvector.segment(j, block_sizes ? (*block_sizes)[block_index] : 1);
                const Scalar block_contribution = eigenvector_block.norm();
                if (block_contribution > contribution_threshold) {
                    std::stringstream internal_ss;
                    internal_ss << std::fixed << std::setprecision(3) << block_contribution;
                    internal_ss << " "
                                << (block_names && !(*block_names)[block_index].empty()
                                                   ? (*block_names)[block_index]
                                                   : (block_sizes ? "Block[" + std::to_string(block_index) + "](" +
                                                                             std::to_string(
                                                                                     (*block_sizes)[block_index]) +
                                                                             ")"
                                                                  : "Element[" + std::to_string(j) + "]"));
                    if (block_sizes) {
                        internal_ss << " [" << eigenvector_block.transpose() << "]";
                    }
                    contributions.emplace_back(internal_ss.str(), block_contribution);
                }
            }

            // Sort eigenvector contributors by significance
            std::sort(contributions.begin(), contributions.end(),
                    [](const std::pair<std::string, Scalar>& lhs, const std::pair<std::string, Scalar>& rhs) {
                        return lhs.second > rhs.second;
                    });

            // Generate diagnostic information
            std::stringstream internal_ss;
            internal_ss << std::scientific << std::setprecision(3) << std::setw(max_index_width) << i << " | "
                        << std::setw(9)
                        << (eigenvalue < nullspace_threshold ? "Nullspace"
                                                             : (eigenvalue < weak_threshold ? "Weak" : "Strong"))
                        << " | Eigenvalue: " << std::setw(9) << eigenvalue << " | Eigenvector contributors: ";
            const std::size_t length_to_contributors = internal_ss.str().size();
            for (std::size_t j = 0; j < contributions.size(); ++j) {
                internal_ss << (j > 0 ? ",\n" + std::string(length_to_contributors, ' ') : "")
                            << contributions[j].first;
            }
            if (contributions.size() < static_cast<std::size_t>(eigenvector.size())) {
                internal_ss << std::fixed << (contributions.size() == 0 ? " None" : "") << ",\n"
                            << std::string(length_to_contributors, ' ')
                            << std::sqrt(static_cast<Scalar>(1) -
                                         std::accumulate(contributions.cbegin(), contributions.cend(),
                                                 static_cast<Scalar>(0),
                                                 [](const Scalar squared_norm,
                                                         const std::pair<std::string, Scalar>& eigenvector_component) {
                                                     return squared_norm +
                                                            eigenvector_component.second * eigenvector_component.second;
                                                 }))
                            << " Remainder";
            }
            internal_ss << "\n";
            if ((eigenvalue < nullspace_threshold && !exclude_nullspaces) ||
                    (eigenvalue >= nullspace_threshold && eigenvalue < weak_threshold && !exclude_weak) ||
                    (eigenvalue >= weak_threshold && !exclude_strong)) {
                ss << internal_ss.str();
            }
        }
        ss << std::string(header_width, '=');
    }
    return ss.str();
}

}

#endif
