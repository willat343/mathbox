#ifndef MATHBOX_TRAITS_HPP
#define MATHBOX_TRAITS_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <type_traits>

namespace math {

template<class ArithmeticType, typename = void>
class ArithmeticTypeTraits {};

template<class ArithmeticType_>
class ArithmeticTypeTraits<ArithmeticType_, typename std::enable_if_t<std::is_arithmetic_v<ArithmeticType_>>> {
public:
    using ArithmeticType = ArithmeticType_;
    using Scalar = ArithmeticType;

    /**
     * @brief Return zero for this type.
     *
     * @return constexpr ArithmeticType
     */
    static constexpr ArithmeticType zero() {
        return static_cast<ArithmeticType>(0);
    }
};

template<class Derived>
class ArithmeticTypeTraits<Derived, typename std::enable_if_t<std::is_base_of_v<Eigen::MatrixBase<Derived>, Derived>>> {
public:
    using ArithmeticType = Derived;
    using Scalar = Derived::Scalar;

    /**
     * @brief Return zero for this type.
     *
     * @return ArithmeticType
     */
    static constexpr ArithmeticType zero() {
        return ArithmeticType::Zero();
    }
};

/**
 * @brief False version of type trait to check if type is a mathbox arithmetic type (of which std::is_arithmetic types
 * are a subset).
 *
 * @tparam T type
 */
template<typename T, typename = void>
struct is_arithmetic_type : std::false_type {};

/**
 * @brief True version of type trait to check if type is a mathbox arithmetic type (of which std::is_arithmetic types
 * are a subset).
 *
 * @tparam T type
 */
template<typename T>
struct is_arithmetic_type<T, std::enable_if_t<std::is_arithmetic_v<T> || std::is_base_of_v<Eigen::MatrixBase<T>, T>>>
    : std::true_type {};

/**
 * @brief Value for type trait to check if type is a mathbox arithmetic type (of which std::is_arithmetic types are a
 * subset).
 *
 * @tparam T type
 */
template<typename T>
constexpr bool is_arithmetic_type_v = is_arithmetic_type<T>::value;

}

#endif
