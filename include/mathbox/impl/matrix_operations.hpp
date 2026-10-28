#ifndef MATHBOX_IMPL_MATRIX_OPERATIONS_HPP
#define MATHBOX_IMPL_MATRIX_OPERATIONS_HPP

#include <algorithm>
#include <cppbox/exceptions.hpp>
#include <limits>

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
    throw_if(m.rows() != v.rows(), "Number of rows must match.", std::invalid_argument);
    throw_if(v.cols() != 1, "Vector must have 1 column.", std::invalid_argument);
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
    throw_if(boundary == 0, "Reorder boundary cannot be 0.", std::invalid_argument);
    const Eigen::Index size = m.rows();
    throw_if(size != m.cols(), "Matrix must be square.", std::invalid_argument);
    throw_if(boundary >= size, "Reorder boundary outside of m matrix.", std::invalid_argument);
    return (Derived(m.rows(), m.cols()) << m.block(boundary, boundary, size - boundary, size - boundary),
            m.block(boundary, 0, size - boundary, boundary), m.block(0, boundary, boundary, size - boundary),
            m.block(0, 0, boundary, boundary))
            .finished();
}

template<typename DerivedMatrix, typename DerivedVector>
    requires(std::is_same_v<typename DerivedMatrix::Scalar, typename DerivedVector::Scalar>)
typename DerivedMatrix::Scalar schur_complement(const Eigen::MatrixBase<DerivedMatrix>& H,
        const Eigen::MatrixBase<DerivedVector>& b, const int upper_block_size,
        Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, Eigen::Dynamic>& H_p,
        Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, 1>& b_p,
        const typename DerivedMatrix::Scalar damping_factor,
        const typename DerivedMatrix::Scalar symmetry_violation_threshold, const bool jacobi_scaling) {
    assert(H.rows() == H.cols() && H.rows() == b.size());
    assert(upper_block_size > 0 && upper_block_size < b.size());
    const int lower_block_size = b.size() - upper_block_size;

    // Apply damping to H_mm, addressing weakly constrained
    const Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, Eigen::Dynamic> H_mm =
            H.topLeftCorner(upper_block_size, upper_block_size);
    const typename DerivedMatrix::Scalar damping =
            damping_factor * H_mm.trace() / static_cast<typename DerivedMatrix::Scalar>(H_mm.rows());

    // Jacobi scaling vector \f$ D^{-1} \f$ where \f$ D = \text{diag}(H_{mm})^{1/2} \f$. A non-positive diagonal entry
    // is replaced by one, leaving it unscaled, so that a non-positive-definite H_mm still fails in the decomposition
    // below rather than a NaN propagating into it.
    Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, 1> D_inv =
            Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, 1>::Ones(upper_block_size);
    if (jacobi_scaling) {
        const Eigen::Array<typename DerivedMatrix::Scalar, Eigen::Dynamic, 1> H_mm_diagonal = H_mm.diagonal();
        D_inv = (H_mm_diagonal > static_cast<typename DerivedMatrix::Scalar>(0))
                        .select(H_mm_diagonal, static_cast<typename DerivedMatrix::Scalar>(1))
                        .rsqrt()
                        .matrix();
    }

    // Create LLT decomposition of the scaled block, D^{-1} (H_mm + damping I) D^{-1}, noting that scaling the damped
    // block is equivalent to adding damping / H_mm(i, i) to the unit diagonal, so damping keeps its unscaled meaning.
    Eigen::LLT<Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, Eigen::Dynamic>> llt(
            D_inv.asDiagonal() *
            (H_mm + damping * Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, Eigen::Dynamic>::Identity(
                                      upper_block_size, upper_block_size)) *
            D_inv.asDiagonal());
    math::check_computation_info(llt.info());

    // Scale the off-diagonal blocks. They are scaled independently rather than transposing one of them, so that any
    // asymmetry of H is preserved into H_p and remains detectable by the symmetry check below.
    const Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, Eigen::Dynamic> H_km_times_D_inv =
            H.bottomLeftCorner(lower_block_size, upper_block_size) * D_inv.asDiagonal();
    const Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, Eigen::Dynamic> D_inv_times_H_mk =
            D_inv.asDiagonal() * H.topRightCorner(upper_block_size, lower_block_size);

    // Solve for the scaled H_mm inverse terms using LLT
    const Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, Eigen::Dynamic> H_mm_inv_times_H_mk =
            llt.solve(D_inv_times_H_mk);
    const Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, 1> H_mm_inv_times_b_m =
            llt.solve(D_inv.asDiagonal() * b.head(upper_block_size));

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

template<typename DerivedMatrix, typename DerivedVector>
    requires(std::is_same_v<typename DerivedMatrix::Scalar, typename DerivedVector::Scalar>)
Eigen::Index square_root_schur_complement(const Eigen::MatrixBase<DerivedMatrix>& J,
        const Eigen::MatrixBase<DerivedVector>& r, const int left_block_size,
        Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, Eigen::Dynamic>& S_p,
        Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, 1>& r_p,
        Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, 1>& pivots_mm,
        Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, 1>& pivots_p,
        const typename DerivedMatrix::Scalar rank_threshold_factor) {
    assert(J.rows() == r.size());
    assert(left_block_size > 0 && left_block_size < J.cols());
    const Eigen::Index right_block_size = J.cols() - left_block_size;

    // Decompose J_m P_m = Q_m R_mm, whose numerical rank p is the number of pivots above the threshold, relative to the
    // maximum pivot, of the conventional tolerance for the dimensions of J
    Eigen::ColPivHouseholderQR<Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, Eigen::Dynamic>> qr_mm(
            J.leftCols(left_block_size));
    qr_mm.setThreshold(numerical_rank_tolerance(J.rows(), J.cols(), rank_threshold_factor));
    const Eigen::Index rank = qr_mm.rank();
    pivots_mm = qr_mm.matrixR().diagonal().cwiseAbs();

    // Rotate [J_k r] by Q_m^T, applied as Householder reflections without forming Q_m, after which the first p rows
    // are satisfied by dx_m for any dx_k. The remaining (marginal) rows give the marginal 0.5 || J~_k dx_k + r~ ||^2,
    // their J_m entries being numerically zero (at most the first pivot below the threshold). The reflectors after the
    // first p only mix the marginal rows, which leaves the marginal unchanged, so they are skipped.
    const Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, Eigen::Dynamic> J_k_r =
            qr_mm.householderQ().setLength(rank).adjoint() *
            (Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, Eigen::Dynamic>(J.rows(),
                     right_block_size + 1)
                            << J.rightCols(right_block_size),
                    r)
                    .finished();
    const Eigen::Index marginal_rows = J.rows() - rank;

    // Decompose J~_k P_p = Q_p R_p, giving S_p = R_p P_p^T and r_p as the leading rows of Q_p^T r~. The matrix R
    // returned by Eigen also holds Householder vectors below its diagonal, so only its upper triangle is taken.
    Eigen::ColPivHouseholderQR<Eigen::Matrix<typename DerivedMatrix::Scalar, Eigen::Dynamic, Eigen::Dynamic>> qr_p(
            J_k_r.bottomLeftCorner(marginal_rows, right_block_size));
    const Eigen::Index S_p_rows = std::min(marginal_rows, right_block_size);
    S_p = qr_p.matrixR().topRows(S_p_rows).template triangularView<Eigen::Upper>();
    S_p = S_p * qr_p.colsPermutation().transpose();
    r_p = (qr_p.householderQ().adjoint() * J_k_r.col(right_block_size).tail(marginal_rows)).head(S_p_rows);
    pivots_p = qr_p.matrixR().diagonal().cwiseAbs();

    // Return the numerical rank of J_m
    return rank;
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

template<typename Derived>
    requires(Derived::RowsAtCompileTime >= 2 && Derived::ColsAtCompileTime == 1)
Eigen::Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, Derived::RowsAtCompileTime - 1>
orthonormal_complement(const Eigen::MatrixBase<Derived>& unit_vector) {
    using Scalar = typename Derived::Scalar;
    constexpr int N = Derived::RowsAtCompileTime;
    throw_if(!unit_vector.isUnitary(), "The vector must be normalised.", std::invalid_argument);
    // k = argmax_i |u_i|
    Eigen::Index k{0};
    unit_vector.cwiseAbs().maxCoeff(&k);
    // w = u + sgn(u_k) e_k
    Eigen::Vector<Scalar, N> w = unit_vector;
    w[k] += unit_vector[k] < static_cast<Scalar>(0) ? static_cast<Scalar>(-1) : static_cast<Scalar>(1);
    // H = I - 2 w w^T / (w^T w)
    const Eigen::Matrix<Scalar, N, N> H =
            Eigen::Matrix<Scalar, N, N>::Identity() - (static_cast<Scalar>(2) / w.squaredNorm()) * w * w.transpose();
    // B = H without its k-th column
    return (Eigen::Matrix<Scalar, N, N - 1>() << H.leftCols(k), H.rightCols(N - 1 - k)).finished();
}

}

#if !MATHBOX_HEADER_ONLY
namespace math {

extern template Eigen::Matrix3d skew_symmetric_cross(const Eigen::MatrixBase<Eigen::Vector3d>& v);

extern template Eigen::Matrix<double, 3, 2> orthonormal_complement(
        const Eigen::MatrixBase<Eigen::Vector3d>& unit_vector);

}
#endif

#endif
