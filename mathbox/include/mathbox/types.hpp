#ifndef MATHBOX_TYPES_HPP
#define MATHBOX_TYPES_HPP

#include <Eigen/Core>

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

}

#endif
