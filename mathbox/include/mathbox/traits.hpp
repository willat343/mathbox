#ifndef MATHBOX_TRAITS_HPP
#define MATHBOX_TRAITS_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <type_traits>

namespace math {

template<class MathType, typename = void>
class MathTypeTraits {};

template<class MathType_>
class MathTypeTraits<MathType_, typename std::enable_if_t<std::is_floating_point_v<MathType_>>> {
public:
    using MathType = MathType_;
    using Scalar = MathType;

    /**
     * @brief Return zero for this type.
     *
     * @return constexpr MathType
     */
    static constexpr MathType zero() {
        return static_cast<MathType>(0);
    }
};

template<class Derived>
class MathTypeTraits<Derived, typename std::enable_if_t<std::is_base_of_v<Eigen::MatrixBase<Derived>, Derived>>> {
public:
    using MathType = Derived;
    using Scalar = typename MathType::Scalar;

    /**
     * @brief Return zero for this type.
     *
     * @return MathType
     */
    static constexpr MathType zero() {
        return MathType::Zero();
    }
};

/**
 * @brief False version of type trait to check if type is a mathbox arithmetic type (of which std::is_arithmetic types
 * are a subset).
 *
 * @tparam T type
 */
template<typename T, typename = void>
struct is_math_type : std::false_type {};

/**
 * @brief True version of type trait to check if type is a mathbox arithmetic type (of which std::is_arithmetic types
 * are a subset).
 *
 * @tparam T type
 */
template<typename T>
struct is_math_type<T, std::enable_if_t<std::is_floating_point_v<T> || std::is_base_of_v<Eigen::MatrixBase<T>, T>>>
    : std::true_type {};

/**
 * @brief Value for type trait to check if type is a mathbox arithmetic type (of which std::is_arithmetic types are a
 * subset).
 *
 * @tparam T type
 */
template<typename T>
constexpr bool is_math_type_v = is_math_type<T>::value;

}

#endif
