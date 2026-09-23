#ifndef MATHBOX_IMPL_MATRIX_OPERATIONS_HPP
#define MATHBOX_IMPL_MATRIX_OPERATIONS_HPP

#include <cppbox/exceptions.hpp>

#include "mathbox/matrix_diagnostics.hpp"
#include "mathbox/matrix_operations.hpp"
#include "mathbox/matrix_properties.hpp"

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

template<typename DerivedMatrix, typename DerivedVector>
    requires(std::is_same_v<typename DerivedMatrix::Scalar, typename DerivedVector::Scalar>)
Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, DerivedMatrix::ColsAtCompileTime>
remove_rows_by_threshold(const Eigen::MatrixBase<DerivedMatrix>& m, const Eigen::MatrixBase<DerivedVector>& v,
        const typename DerivedMatrix::Scalar threshold) {
    throw_if(m.rows() != v.rows(), "Number of rows must match.");
    throw_if(v.cols() != 1, "Vector must have 1 column.");
    Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, DerivedMatrix::ColsAtCompileTime> m_reduced(m.rows(),
            m.cols());
    int m_reduced_rows{0};
    for (int r = 0; r < m.rows(); ++r) {
        if (v[r] >= threshold) {
            m_reduced.row(m_reduced_rows) = m.row(r);
            ++m_reduced_rows;
        }
    }
    m_reduced.conservativeResize(m_reduced_rows, m.cols());
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
        const double symmetry_violation_threshold, const bool jacobi_scaling) {
    assert(H.rows() == H.cols() && H.rows() == b.size());
    assert(upper_block_size > 0 && upper_block_size < b.size());
    const int lower_block_size = b.size() - upper_block_size;

    // Apply damping to H_mm, addressing weakly constrained
    const Eigen::MatrixXd H_mm = H.topLeftCorner(upper_block_size, upper_block_size);
    const double damping = damping_factor * H_mm.trace() / static_cast<double>(H_mm.rows());

    // Jacobi scaling vector \f$ D^{-1} \f$ where \f$ D = \text{diag}(H_{mm})^{1/2} \f$. A non-positive diagonal entry
    // is replaced by one, leaving it unscaled, so that a non-positive-definite H_mm still fails in the decomposition
    // below rather than a NaN propagating into it.
    Eigen::VectorXd D_inv = Eigen::VectorXd::Ones(upper_block_size);
    if (jacobi_scaling) {
        const Eigen::ArrayXd H_mm_diagonal = H_mm.diagonal();
        D_inv = (H_mm_diagonal > 0.0).select(H_mm_diagonal, 1.0).rsqrt().matrix();
    }

    // Create LLT decomposition of the scaled block, D^{-1} (H_mm + damping I) D^{-1}, noting that scaling the damped
    // block is equivalent to adding damping / H_mm(i, i) to the unit diagonal, so damping keeps its unscaled meaning.
    Eigen::LLT<Eigen::MatrixXd> llt(D_inv.asDiagonal() *
                                    (H_mm + damping * Eigen::MatrixXd::Identity(upper_block_size, upper_block_size)) *
                                    D_inv.asDiagonal());
    math::check_computation_info(llt.info());

    // Scale the off-diagonal blocks. They are scaled independently rather than transposing one of them, so that any
    // asymmetry of H is preserved into H_p and remains detectable by the symmetry check below.
    const Eigen::MatrixXd H_km_times_D_inv =
            H.bottomLeftCorner(lower_block_size, upper_block_size) * D_inv.asDiagonal();
    const Eigen::MatrixXd D_inv_times_H_mk = D_inv.asDiagonal() * H.topRightCorner(upper_block_size, lower_block_size);

    // Solve for the scaled H_mm inverse terms using LLT
    const Eigen::MatrixXd H_mm_inv_times_H_mk = llt.solve(D_inv_times_H_mk);
    const Eigen::MatrixXd H_mm_inv_times_b_m = llt.solve(D_inv.asDiagonal() * b.head(upper_block_size));

    // Compute H_p and b_p
    H_p = H.bottomRightCorner(lower_block_size, lower_block_size) - H_km_times_D_inv * H_mm_inv_times_H_mk;
    b_p = b.tail(lower_block_size) - H_km_times_D_inv * H_mm_inv_times_b_m;

    // H may not be exactly symmetric due to numerical precision, so enforce symmetry
    throw_if(math::relative_asymmetry(H_p) >= symmetry_violation_threshold, "Symmetry threshold violated for H_p.");
    H_p = math::make_symmetric(H_p);
    assert(math::relative_asymmetry(H_p) < symmetry_violation_threshold);

    // Return the damping
    return damping;
}

template<typename Derived>
    requires(Derived::RowsAtCompileTime == 3 && Derived::ColsAtCompileTime == 1)
constexpr inline Eigen::Matrix<typename Derived::Scalar, 3, 3> skew_symmetric_cross(
        const Eigen::MatrixBase<Derived>& v) {
    return (Eigen::Matrix<typename Derived::Scalar, 3, 3>() << static_cast<typename Derived::Scalar>(0), -v[2], v[1],
            v[2], static_cast<typename Derived::Scalar>(0), -v[0], -v[1], v[0],
            static_cast<typename Derived::Scalar>(0))
            .finished();
}

}

#if !MATHBOX_HEADER_ONLY
namespace math {

extern template Eigen::Matrix3d skew_symmetric_cross(const Eigen::MatrixBase<Eigen::Vector3d>& v);

}
#endif

#endif
