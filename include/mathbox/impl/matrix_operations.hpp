#ifndef MATHBOX_IMPL_MATRIX_OPERATIONS_HPP
#define MATHBOX_IMPL_MATRIX_OPERATIONS_HPP

#include <cppbox/exceptions.hpp>

#include "mathbox/matrix_diagnostics.hpp"
#include "mathbox/matrix_operations.hpp"

namespace math {

template<typename Derived>
constexpr Derived cumulative_row_left_sum(const Eigen::DenseBase<Derived>& m) {
    if (m.cols() == 0) {
        return m.derived();
    }
    Derived m_cumulative(m.rows(), m.cols());
    m_cumulative.col(m_cumulative.cols() - 1) = m.col(m_cumulative.cols() - 1);
    for (int c = m_cumulative.cols() - 2; c >= 0; --c) {
        m_cumulative.col(c) = m_cumulative.col(c + 1) + m.col(c);
    }
    return m_cumulative;
}

template<typename Derived>
constexpr Derived cumulative_col_top_sum(const Eigen::DenseBase<Derived>& m) {
    if (m.rows() == 0) {
        return m.derived();
    }
    Derived m_cumulative(m.rows(), m.cols());
    m_cumulative.row(m_cumulative.rows() - 1) = m.row(m_cumulative.rows() - 1);
    for (int r = m_cumulative.rows() - 2; r >= 0; --r) {
        m_cumulative.row(r) = m_cumulative.row(r + 1) + m.row(r);
    }
    return m_cumulative;
}

template<typename Derived>
constexpr typename Derived::PlainObject make_symmetric(const Eigen::MatrixBase<Derived>& m) {
    return 0.5 * (m + m.transpose());
}

template<typename Derived>
void make_symmetric_inplace(Eigen::MatrixBase<Derived>& m) {
    m = 0.5 * (m + m.transpose());
}

template<typename DerivedMatrix, typename DerivedVector>
    requires(std::is_same_v<typename DerivedMatrix::Scalar, typename DerivedVector::Scalar>)
Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, DerivedMatrix::ColsAtCompileTime>
remove_rows_by_threshold(const Eigen::MatrixBase<DerivedMatrix>& m, const Eigen::MatrixBase<DerivedVector>& v,
        const typename DerivedMatrix::Scalar threshold) {
    throw_if(m.rows() != v.rows(), "Number of rows must match.");
    throw_if(v.cols() != 1, "Vector must have 1 column.");
    Eigen::MatrixXd m_reduced(m.rows(), m.cols());
    int m_reduced_rows{0};
    for (int r = 0; r < m.rows(); ++r) {
        if (v[r] >= threshold) {
            m_reduced.row(m_reduced_rows) = m.row(r);
            ++m_reduced_rows;
        }
    }
    m_reduced.resize(m_reduced_rows, m.cols());
    return m_reduced;
}

template<typename Derived>
Derived reorder_symmetric_matrix(const Eigen::MatrixBase<Derived>& m, const Eigen::Index boundary) {
    throw_if(boundary == 0, "Reorder boundary cannot be 0.");
    const Eigen::Index size = m.rows();
    throw_if(size != m.cols(), "Matrix must be square.");
    throw_if(boundary >= size, "Reorder boundary outside of m matrix.");
    return (Derived(m.rows(), m.cols()) << m.block(boundary, boundary, size - boundary, size - boundary),
            m.block(boundary, 0, size - boundary, boundary), m.block(0, boundary, boundary, size - boundary),
            m.block(0, 0, boundary, boundary))
            .finished();
}

inline double schur_complement(const Eigen::Ref<const Eigen::MatrixXd>& H, const Eigen::Ref<const Eigen::VectorXd> b,
        const int upper_block_size, Eigen::MatrixXd& H_p, Eigen::VectorXd& b_p, const double damping_factor,
        const double symmetry_violation_threshold) {
    assert(H.rows() == H.cols() && H.rows() == b.size());
    assert(upper_block_size > 0 && upper_block_size < b.size() - 1);
    const int lower_block_size = b.size() - upper_block_size;

    // Apply damping to H_mm, addressing weakly constrained
    const Eigen::MatrixXd H_mm = H.topLeftCorner(upper_block_size, upper_block_size);
    const double damping = damping_factor * H_mm.trace() / static_cast<double>(H_mm.rows());

    // Create LLT decomposition
    Eigen::LLT<Eigen::MatrixXd> llt(H_mm + damping * Eigen::MatrixXd::Identity(upper_block_size, upper_block_size));
    math::check_computation_info(llt.info());

    // Solve for the H_mm inverse terms using LLT
    const Eigen::MatrixXd H_mm_inv_times_H_mk = llt.solve(H.topRightCorner(upper_block_size, lower_block_size));
    const Eigen::MatrixXd H_mm_inv_times_b_m = llt.solve(b.head(upper_block_size));

    // Compute H_p and b_p
    H_p = H.bottomRightCorner(lower_block_size, lower_block_size) -
          H.bottomLeftCorner(lower_block_size, upper_block_size) * H_mm_inv_times_H_mk;
    b_p = b.tail(lower_block_size) - H.bottomLeftCorner(lower_block_size, upper_block_size) * H_mm_inv_times_b_m;

    // H may not be exactly symmetric due to numerical precision, so enforce symmetry
    throw_if((H_p - H_p.transpose()).norm() / H_p.norm() >= symmetry_violation_threshold,
            "Symmetry threshold violated for H_p.");
    math::make_symmetric_inplace(H_p);
    assert((H_p - H_p.transpose()).norm() / H_p.norm() < symmetry_violation_threshold);

    // Return the damping
    return damping;
}

template<typename Scalar>
inline Eigen::Matrix<Scalar, 3, 3> skew_symmetric_cross(const Eigen::Matrix<Scalar, 3, 1>& v) {
    return (Eigen::Matrix<Scalar, 3, 3>() << static_cast<Scalar>(0), -v[2], v[1], v[2], static_cast<Scalar>(0), -v[0],
            -v[1], v[0], static_cast<Scalar>(0))
            .finished();
}

}

#if !MATHBOX_HEADER_ONLY
namespace math {

extern template Eigen::Matrix3d skew_symmetric_cross(const Eigen::Vector3d& v);

}
#endif

#endif
