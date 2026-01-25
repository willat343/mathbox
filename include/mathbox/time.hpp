#ifndef MATHBOX_TIME_HPP
#define MATHBOX_TIME_HPP

#include <Eigen/Core>
#include <chrono>
#include <cppbox/time.hpp>
#include <vector>

namespace math {

/**
 * @brief Generate a linearly spaced dynamically-sized matrix from an interpolation step. The number of points are
 * chosen such that the actual interpolation step is greater than or equal to `step`.
 *
 * The first element will be `start` and the last element will be `end`.
 *
 * Overload of math::lin_spaced_vector for Time and Duration.
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
