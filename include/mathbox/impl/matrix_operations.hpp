#ifndef MATHBOX_IMPL_MATRIX_OPERATIONS_HPP
#define MATHBOX_IMPL_MATRIX_OPERATIONS_HPP

#include "mathbox/matrix_operations.hpp"

namespace math {

template<typename Derived>
constexpr Derived cumulative_row_left_sum(const Eigen::DenseBase<Derived>& m) {
    if (m.cols() == 0) {
        return m.derived();
    }
    Derived m_cumulative;
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
    Derived m_cumulative;
    m_cumulative.row(m_cumulative.rows() - 1) = m.row(m_cumulative.rows() - 1);
    for (int r = m_cumulative.rows() - 2; r >= 0; --r) {
        m_cumulative.row(r) = m_cumulative.row(r + 1) + m.row(r);
    }
    return m_cumulative;
}

}

#endif
