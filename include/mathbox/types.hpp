#ifndef MATHBOX_TYPES_HPP
#define MATHBOX_TYPES_HPP

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <type_traits>

#include "mathbox/traits.hpp"

namespace math {

#define MATHBOX_MAP_TYPES(Name)         \
    using Name##Map = Eigen::Map<Name>; \
    using Name##ConstMap = Eigen::Map<const Name>;

#define MATHBOX_REF_TYPES(Name)         \
    using Name##Ref = Eigen::Map<Name>; \
    using Name##ConstRef = const Eigen::Ref<const Name>&;

#define MATHBOX_MATRIX_TYPES(Name) \
    MATHBOX_MAP_TYPES(Name)        \
    MATHBOX_REF_TYPES(Name)

/**
 * @brief Orientation type.
 *
 * @tparam D
 */
template<int D>
    requires(is_2d_or_3d<D>)
using Orientation = std::conditional_t<D == 2, Eigen::Rotation2D<double>, Eigen::Quaternion<double>>;

/**
 * @brief Pose type.
 *
 * @tparam D dimension
 */
template<int D>
    requires(is_2d_or_3d<D>)
using Pose = Eigen::Transform<double, D, Eigen::Isometry>;

/**
 * @brief Position type.
 *
 * @tparam D
 */
template<int D>
    requires(is_2d_or_3d<D>)
using Position = Eigen::Vector<double, D>;

}

#endif
