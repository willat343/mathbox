#ifndef MATHBOX_TRAITS_HPP
#define MATHBOX_TRAITS_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cppbox/time.hpp>
#include <type_traits>

namespace math {

/**
 * @brief Trait for checking the dimension D is 2D or 3D, i.e. the integer is 2 or 3.
 *
 * Example usage:
 * ```
 * template<int D>
 *     requires is_2d_or_3d<D>
 * void foo();
 * ```
 *
 * @tparam D
 */
template<int D>
constexpr bool is_2d_or_3d = (D == 2 || D == 3);

template<typename T>
concept IsMatrix = requires(T t) {
    { static_cast<const Eigen::MatrixBase<T>&>(t) };
};

template<typename T, typename Scalar>
concept IsMatrixOf = IsMatrix<T> && std::is_same_v<typename T::Scalar, Scalar>;

template<typename T>
concept IsMatrixf = IsMatrixOf<T, float>;

template<typename T>
concept IsMatrixd = IsMatrixOf<T, double>;

template<typename T>
concept IsSparseMatrix = requires(T t) {
    { static_cast<const Eigen::SparseMatrixBase<T>&>(t) };
};

template<typename T, typename Scalar>
concept IsSparseMatrixOf = IsSparseMatrix<T> && std::is_same_v<typename T::Scalar, Scalar>;

template<typename T>
concept IsSparseMatrixf = IsSparseMatrixOf<T, float>;

template<typename T>
concept IsSparseMatrixd = IsSparseMatrixOf<T, double>;

template<typename T>
concept IsRowMajorMatrix = IsMatrix<T> && T::IsRowMajor == 1;

template<typename T, typename Scalar>
concept IsRowMajorMatrixOf = IsRowMajorMatrix<T> && std::is_same_v<typename T::Scalar, Scalar>;

template<typename T>
concept IsRowMajorMatrixf = IsRowMajorMatrixOf<T, float>;

template<typename T>
concept IsRowMajorMatrixd = IsRowMajorMatrixOf<T, double>;

template<typename T>
concept IsRowMajorSparseMatrix = IsSparseMatrix<T> && T::IsRowMajor == 1;

template<typename T, typename Scalar>
concept IsRowMajorSparseMatrixOf = IsRowMajorSparseMatrix<T> && std::is_same_v<typename T::Scalar, Scalar>;

template<typename T>
concept IsRowMajorSparseMatrixf = IsRowMajorSparseMatrixOf<T, float>;

template<typename T>
concept IsRowMajorSparseMatrixd = IsRowMajorSparseMatrixOf<T, double>;

template<typename T>
concept IsVector = IsMatrix<T> && T::IsVectorAtCompileTime == 1;

template<typename T, typename Scalar>
concept IsVectorOf = IsVector<T> && std::is_same_v<typename T::Scalar, Scalar>;

template<typename T>
concept IsVectorf = IsVectorOf<T, float>;

template<typename T>
concept IsVectord = IsVectorOf<T, double>;

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
