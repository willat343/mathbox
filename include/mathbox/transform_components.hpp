#ifndef MATHBOX_TRANSFORM_COMPONENTS_HPP
#define MATHBOX_TRANSFORM_COMPONENTS_HPP

#include <Eigen/Core>

namespace math {

/**
 * @brief Precompute and store the components of a transformation. The purpose of this class is primarily to avoid
 * repeated computation, when these components are needed many times.
 * 
 * @tparam Scalar 
 */
template<typename Scalar>
class TransformComponents {
public:
    explicit TransformComponents(const Eigen::Matrix<Scalar, 4, 4>& transform);

    explicit TransformComponents(const Eigen::Matrix<Scalar, 3, 1>& r, const Eigen::Matrix<Scalar, 3, 1>& t);

    inline const Eigen::Matrix<Scalar, 4, 4>& T() const;
    inline const Eigen::Matrix<Scalar, 3, 3>& R() const;
    inline const Eigen::Matrix<Scalar, 3, 1>& r() const;
    inline const Eigen::Matrix<Scalar, 3, 1>& t() const;
    inline const Scalar a() const;
    inline const Scalar sina() const;
    inline const Scalar cosa() const;

private:
    Eigen::Matrix<Scalar, 4, 4> T_;
    Eigen::Matrix<Scalar, 3, 3> R_;
    Eigen::Matrix<Scalar, 3, 1> r_;
    Eigen::Matrix<Scalar, 3, 1> t_;
    Scalar a_;
    Scalar sina_;
    Scalar cosa_;
};

}

#include "mathbox/transform_components.hpp"

#endif
