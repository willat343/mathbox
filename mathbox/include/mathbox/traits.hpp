#ifndef MATHBOX_TRAITS_HPP
#define MATHBOX_TRAITS_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <cppbox/time.hpp>
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
        static_assert(MathType::SizeAtCompileTime != Eigen::Dynamic,
                "MathTypeTraits::zero() cannot be called for dynamic-sized matrices.");
        return MathType::Zero();
    }

    static constexpr MathType zero(const int rows, const int cols) {
        return MathType::Zero(rows, cols);
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

/**
 * @brief Math type concept.
 *
 * @tparam T
 */
template<typename T>
concept IsMathType = is_math_type_v<T>;

template<typename T, typename = void>
struct is_independent_variable_type : std::false_type {};

template<typename T>
struct is_independent_variable_type<T, std::enable_if_t<std::is_floating_point_v<T> || cppbox::is_time_point_v<T>>>
    : std::true_type {};

template<typename T>
constexpr bool is_independent_variable_type_v = is_independent_variable_type<T>::value;

template<typename T>
concept IsIndependentVariableType = is_independent_variable_type_v<T>;

}

#endif
