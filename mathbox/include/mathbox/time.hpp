#ifndef MATHBOX_TIME_HPP
#define MATHBOX_TIME_HPP

#include <chrono>

namespace math {

/**
 * @brief False version of type trait to check if type is a `std::chrono::Duration<Rep, Period>`.
 *
 * @tparam T type
 */
template<typename T>
struct is_duration : std::false_type {};

/**
 * @brief True version of type trait to check if type is a `std::chrono::Duration<Rep, Period>`.
 *
 * @tparam Rep
 * @tparam Period
 */
template<class Rep, class Period>
struct is_duration<std::chrono::duration<Rep, Period>> : std::true_type {};

/**
 * @brief Value for type trait to check if type is a `std::chrono::Duration<Rep, Period>`.
 *
 * @tparam T type
 */
template<typename T>
constexpr bool is_duration_v = is_duration<T>::value;

/**
 * @brief False version of type trait to check if type is a `std::chrono::time_point<Clock, Duration>`.
 *
 * @tparam T type
 */
template<typename T>
struct is_time_point : std::false_type {};

/**
 * @brief True version of type trait to check if type is a `std::chrono::time_point<Clock, Duration>`.
 *
 * @tparam Clock
 * @tparam Duration
 */
template<class Clock, class Duration>
struct is_time_point<std::chrono::time_point<Clock, Duration>> : std::true_type {};

/**
 * @brief Value for type trait to check if type is a `std::chrono::time_point<Clock, Duration>`.
 *
 * @tparam T type
 */
template<typename T>
constexpr bool is_time_point_v = is_time_point<T>::value;

/**
 * @brief Compute the fraction (ratio) between two durations as a Scalar value (avoiding possible integer division).
 *
 * @tparam Scalar
 * @tparam Duration
 * @param numerator
 * @param denominator
 * @return Scalar
 */
template<typename Scalar = double, class Duration = Scalar>
constexpr Scalar fraction(const Duration numerator, const Duration denominator);

/**
 * @brief Convert seconds to duration.
 *
 * @tparam Duration type satisfying `std::is_arithemetic_v` or `is_duration_v`
 * @tparam Scalar
 * @param seconds
 * @return Duration
 */
template<class Duration, typename Scalar = double>
constexpr Duration to_duration(const Scalar seconds);

/**
 * @brief Convert time (since clock epoch) or duration to scalar seconds.
 *
 * @tparam Scalar default double allows use of function without specifying template parameters
 * @tparam TimeOrDuration type satisfying `std::is_arithemetic_v`, `is_time_point_v` or `is_duration_v`
 * @param time_or_duration
 * @return Scalar
 */
template<typename Scalar = double, class TimeOrDuration = Scalar>
constexpr Scalar to_sec(const TimeOrDuration& time_or_duration);

/**
 * @brief Convert seconds to time (since clock epoch).
 *
 * @tparam Time type satisfying `std::is_arithemetic_v` or `is_time_point_v`
 * @tparam Scalar
 * @param seconds
 * @return Time
 */
template<class Time, typename Scalar = double>
constexpr Time to_time(const Scalar seconds);

}

#include "mathbox/impl/time.hpp"

#endif
