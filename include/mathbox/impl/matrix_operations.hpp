#ifndef MATHBOX_IMPL_MATRIX_OPERATIONS_HPP
#define MATHBOX_IMPL_MATRIX_OPERATIONS_HPP

#include <cppbox/exceptions.hpp>

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
