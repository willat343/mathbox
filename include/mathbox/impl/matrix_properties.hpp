#ifndef MATHBOX_IMPL_MATRIX_PROPERTIES_HPP
#define MATHBOX_IMPL_MATRIX_PROPERTIES_HPP

#include "mathbox/matrix_properties.hpp"

namespace math {

template<typename Derived>
inline bool is_skew_symmetric(const Eigen::DenseBase<Derived>& m) {
    return m.isApprox(-m.transpose());
}

}

#endif
