#ifndef MATHBOX_IMPL_TRANSFORM_COMPONENTS_HPP
#define MATHBOX_IMPL_TRANSFORM_COMPONENTS_HPP

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "mathbox/transform_components.hpp"

namespace math {

template<typename Scalar>
TransformComponents<Scalar>::TransformComponents(const Eigen::Matrix<Scalar, 4, 4>& transform)
    : T_(transform) {
    R_ = T_.template block<3, 3>(0, 0);
    const Eigen::AngleAxis<Scalar> angleaxis{R_};
    a_ = angleaxis.angle();
    sina_ = std::sin(a_);
    cosa_ = std::cos(a_);
    r_ = angleaxis.axis() * a_;
    t_ = T_.template block<3, 1>(0, 3);
}

template<typename Scalar>
TransformComponents<Scalar>::TransformComponents(const Eigen::Matrix<Scalar, 3, 1>& r,
        const Eigen::Matrix<Scalar, 3, 1>& t)
    : r_(r),
      t_(t) {
    a_ = r_.norm();
    const Eigen::AngleAxis<Scalar> angleaxis{a_, r_ / a_};
    R_ = angleaxis.toRotationMatrix();
    sina_ = std::sin(a_);
    cosa_ = std::cos(a_);
    T_ << R_, t_, Eigen::Matrix<Scalar, 1, 3>::Zero(), static_cast<Scalar>(1);
}

template<typename Scalar>
inline const Eigen::Matrix<Scalar, 4, 4>& TransformComponents<Scalar>::T() const {
    return T_;
}

template<typename Scalar>
inline const Eigen::Matrix<Scalar, 3, 3>& TransformComponents<Scalar>::R() const {
    return R_;
}

template<typename Scalar>
inline const Eigen::Matrix<Scalar, 3, 1>& TransformComponents<Scalar>::r() const {
    return r_;
}

template<typename Scalar>
inline const Eigen::Matrix<Scalar, 3, 1>& TransformComponents<Scalar>::t() const {
    return t_;
}

template<typename Scalar>
inline const Scalar TransformComponents<Scalar>::a() const {
    return a_;
}

template<typename Scalar>
inline const Scalar TransformComponents<Scalar>::sina() const {
    return sina_;
}

template<typename Scalar>
inline const Scalar TransformComponents<Scalar>::cosa() const {
    return cosa_;
}

}

#if !MATHBOX_HEADER_ONLY
namespace math {

extern template class TransformComponents<double>;

}
#endif

#endif
