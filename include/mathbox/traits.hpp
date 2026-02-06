#ifndef MATHBOX_TRAITS_HPP
#define MATHBOX_TRAITS_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cppbox/time.hpp>
#include <cppbox/traits.hpp>
#include <type_traits>

namespace cppbox {

template<typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
struct const_ref<Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>> {
    using type = const Eigen::Ref<const Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>>&;
};

template<typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
struct ref<Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>> {
    using type = Eigen::Ref<Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>>;
};

template<typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
struct ref<const Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>> {
    using type = Eigen::Ref<const Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>>;
};

template<typename PlainObjectType>
struct remove_ref<Eigen::Ref<PlainObjectType>> {
    using type = PlainObjectType;
};

}

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

template<typename T, typename Enable = void>
struct remove_map {
    using type = T;
};

template<typename T, typename Enable = void>
using remove_map_t = remove_map<T, Enable>::type;

template<typename PlainObjectType>
struct remove_map<Eigen::Map<PlainObjectType>> {
    using type = PlainObjectType;
};

template<typename T, typename Enable = void>
struct remove_ref_or_map {
    using type = std::remove_reference_t<T>;
};

template<typename T, typename Enable = void>
using remove_ref_or_map_t = remove_ref_or_map<T, Enable>::type;

template<typename PlainObjectType>
struct remove_ref_or_map<Eigen::Ref<PlainObjectType>> {
    using type = PlainObjectType;
};

template<typename PlainObjectType>
struct remove_ref_or_map<Eigen::Map<PlainObjectType>> {
    using type = PlainObjectType;
};

template<typename T, typename Base, typename Enable = void>
struct is_ref_or_map : std::false_type {};

template<typename T, typename Base, typename Enable = void>
constexpr inline bool is_ref_or_map_v = is_ref_or_map<T, Base, Enable>::value;

template<typename T, typename Base>
struct is_ref_or_map<T, Base,
        std::enable_if_t<std::is_same_v<T, cppbox::ref_t<Base>> || std::is_same_v<T, Eigen::Map<Base>>>>
    : std::true_type {};

/**
 * @brief Reference or Eigen non-const map (Eigen::Map<...>) concept.
 *
 * Types that obey this constraint are typically passed with perfect forwarding, e.g.
 * ```
 * template<IsRefOrMap<Base> T>
 * void foo(T&& lie_group);
 * ```
 *
 * @tparam T type
 * @tparam Base base type
 */
template<typename T, typename Base>
concept IsRefOrMap = is_ref_or_map_v<T, Base>;

template<typename T, typename Base, typename Enable = void>
struct is_same_or_const_map : std::false_type {};

template<typename T, typename Base, typename Enable = void>
constexpr inline bool is_same_or_const_map_v = is_same_or_const_map<T, Base, Enable>::value;

template<typename T, typename Base>
struct is_same_or_const_map<T, Base,
        std::enable_if_t<std::is_same_v<T, Base> || std::is_same_v<T, Eigen::Map<const Base>>>> : std::true_type {};

/**
 * @brief Same or Eigen const map (Eigen::Map<const ...>) concept.
 *
 * Types that obey this constraint are typically passed as const lvalue references, e.g.
 * ```
 * template<IsSameOrConstMap<Base> T>
 * void foo(const T& lie_group);
 * ```
 *
 * This concept is stricter than `IsSameOrAnyMap` since it does not allow Eigen non-const maps (Eigen::Map<...>) even
 * though typical usages passes these types by const reference. Use `IsSameOrAnyMap` for more flexibility.
 *
 * @tparam T type
 * @tparam Base base type
 */
template<typename T, typename Base>
concept IsSameOrConstMap = is_same_or_const_map_v<T, Base>;

template<typename T, typename Base, typename Enable = void>
struct is_same_or_any_map : std::false_type {};

template<typename T, typename Base, typename Enable = void>
constexpr inline bool is_same_or_any_map_v = is_same_or_any_map<T, Base, Enable>::value;

template<typename T, typename Base>
struct is_same_or_any_map<T, Base,
        std::enable_if_t<std::is_same_v<T, Base> || std::is_same_v<T, Eigen::Map<Base>> ||
                         std::is_same_v<T, Eigen::Map<const Base>>>> : std::true_type {};

/**
 * @brief Same or any Eigen non-const map (Eigen::Map<...> or Eigen::Map<const ...>) concept. By collecting these
 * types, this concept allows for more flexible passing of parameters that are meant to be passed by const reference.
 *
 * Types that obey this constraint are typically passed as const lvalue references, e.g.
 * ```
 * template<IsSameOrAnyMap<Base> T>
 * void foo(const T& lie_group);
 * ```
 *
 * Note, however, that in the case of Eigen non-const maps, the constness could in theory be cast away and the object
 * modified.
 *
 * @tparam T type
 * @tparam Base base type
 */
template<typename T, typename Base>
concept IsSameOrAnyMap = is_same_or_any_map_v<T, Base>;

}

#endif
