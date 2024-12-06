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
 * @brief Convert seconds to duration.
 *
 * @tparam Duration type satisfying `std::is_arithemetic_v` or `is_duration_v`
 * @tparam Scalar
 * @param seconds
 * @return Duration
 */
template<class Duration, typename Scalar = double>
Duration to_duration(const Scalar seconds);

/**
 * @brief Convert time (since clock epoch) or duration to scalar seconds.
 *
 * @tparam TimeOrDuration type satisfying `std::is_arithemetic_v`, `is_time_point_v` or `is_duration_v`
 * @tparam Scalar
 * @param time_or_duration
 * @return Scalar
 */
template<typename Scalar, class TimeOrDuration>
Scalar to_sec(const TimeOrDuration& time_or_duration);

/**
 * @brief Convert time (since clock epoch) or duration to double seconds. This function overload allows the user to
 * avoid explicitly stating the Scalar output type as double.
 *
 * @tparam TimeOrDuration type satisfying `std::is_arithemetic_v`, `is_time_point_v` or `is_duration_v`
 * @param time_or_duration
 * @return double
 */
template<class TimeOrDuration>
double to_sec(const TimeOrDuration& time_or_duration);

/**
 * @brief Convert seconds to time (since clock epoch).
 *
 * @tparam Time type satisfying `std::is_arithemetic_v` or `is_time_point_v`
 * @tparam Scalar
 * @param seconds
 * @return Time
 */
template<class Time, typename Scalar = double>
Time to_time(const Scalar seconds);

}

#include "mathbox/impl/time.hpp"

#endif
