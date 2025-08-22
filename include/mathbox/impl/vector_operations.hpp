#ifndef MATHBOX_IMPL_VECTOR_OPERATIONS_HPP
#define MATHBOX_IMPL_VECTOR_OPERATIONS_HPP

#include "mathbox/vector_operations.hpp"

namespace math {

template<typename Scalar>
Eigen::Matrix<Scalar, Eigen::Dynamic, 1> lin_spaced(const Scalar step, const Scalar start, const Scalar end) {
    const Eigen::Index size = static_cast<Eigen::Index>((end - start) / step) + 1;
    return Eigen::Matrix<Scalar, Eigen::Dynamic, 1>::LinSpaced(size, start, end);
}

template<typename Scalar>
Eigen::Matrix<Scalar, Eigen::Dynamic, 1> range(const Scalar step, const Scalar start, const Scalar end) {
    return lin_spaced(step, start, end - std::fmod<Scalar, Scalar>(end - start, step));
}

}

#if !MATHBOX_HEADER_ONLY
namespace math {

extern template Eigen::Matrix<double, Eigen::Dynamic, 1> lin_spaced<double>(const double, const double, const double);

extern template Eigen::Matrix<double, Eigen::Dynamic, 1> range<double>(const double, const double, const double);

}
#endif

#endif
