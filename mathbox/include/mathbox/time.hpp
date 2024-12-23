#ifndef MATHBOX_TIME_HPP
#define MATHBOX_TIME_HPP

#include <Eigen/Core>
#include <chrono>
#include <vector>

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
 * @brief Generate a linearly spaced dynamically-sized matrix from an interpolation step. The number of points are
 * chosen such that the actual interpolation step is greater than or equal to `step`.
 *
 * The first element will be `start` and the last element will be `end`.
 *
 * Overload of math::lin_spaced for Time and Duration.
 *
 * @tparam Time
 * @tparam Time::Duration
 * @param step
 * @param start
 * @param end
 * @return std::vector<Time>
 */
template<class Time, class Duration>
std::vector<Time> lin_spaced(const Duration step, const Time start, const Time end);

/**
 * @brief Generate a range from `start` to `end` incrementing by `step`.
 *
 * The first element will be `start` and the last element will be `end` minus the modulus of (`end` - `start`) and
 * `step`. It is thus at most `end`.
 *
 * Overload of math::range for Time and Duration.
 *
 * @tparam Time
 * @tparam Duration
 * @param step
 * @param start
 * @param end
 * @return std::vector<Time>
 */
template<class Time, class Duration>
std::vector<Time> range(const Duration step, const Time start, const Time end);

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
 * @brief Convert times (since clock epoch) or durations to scalar seconds.
 *
 * @tparam Scalar
 * @tparam TimeOrDuration
 * @param times_or_durations
 * @return std::vector<Scalar>
 */
template<typename Scalar = double, class TimeOrDuration = Scalar>
std::vector<Scalar> to_secs(const std::vector<TimeOrDuration>& times_or_durations);

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

/**
 * @brief Convert vector of seconds to vector of times (since clock epoch).
 *
 * @tparam Time Time type satisfying `std::is_arithemetic_v` or `is_time_point_v`
 * @tparam Scalar
 * @param seconds
 * @return std::vector<Time>
 */
template<class Time, typename Scalar = double>
std::vector<Time> to_times(const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& seconds);

/**
 * @brief Overload of `to_times` for row vectors.
 *
 * @tparam Time
 * @tparam Scalar
 * @param seconds
 * @return std::vector<Time>
 */
template<class Time, typename Scalar = double>
std::vector<Time> to_times(const Eigen::Matrix<Scalar, 1, Eigen::Dynamic>& seconds);

}

#include "mathbox/impl/time.hpp"

#endif
