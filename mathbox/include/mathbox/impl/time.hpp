#ifndef MATHBOX_IMPL_TIME_HPP
#define MATHBOX_IMPL_TIME_HPP

#include <chrono>

#include "mathbox/time.hpp"
#include "mathbox/vector_operations.hpp"

namespace math {

template<class Time, class Duration>
std::vector<Time> lin_spaced(const Duration step, const Time start, const Time end) {
    Eigen::Matrix<double, Eigen::Dynamic, 1> lin_spaced_vector =
            lin_spaced<double>(math::to_sec(step), math::to_sec(start), math::to_sec(end));
    const std::size_t size = static_cast<std::size_t>(lin_spaced_vector.size());
    std::vector<Time> lin_spaced_times(size);
    // Eigen 3.3 does not support iterators so cannot use std::transform.
    for (std::size_t i = 0; i < size; ++i) {
        lin_spaced_times[i] = math::to_time<Time>(lin_spaced_vector[i]);
    }
    return lin_spaced_times;
}

template<class Time, class Duration>
std::vector<Time> range(const Duration step, const Time start, const Time end) {
    const std::size_t size = static_cast<std::size_t>((end - start) / step) + 1;
    std::vector<Time> range_times(size);
    for (std::size_t i = 0; i < size; ++i) {
        range_times[i] = start + i * step;
    }
    return range_times;
}

template<class Duration, typename Scalar>
constexpr inline Duration to_duration(const Scalar seconds) {
    static_assert(std::is_floating_point_v<Duration> || is_duration_v<Duration>,
            "Duration must be of floating point scalar or duration type.");
    if constexpr (std::is_floating_point_v<Duration>) {
        return static_cast<Duration>(seconds);
    } else if constexpr (is_duration_v<Duration>) {
        return std::chrono::duration_cast<Duration>(std::chrono::duration<Scalar>(seconds));
    } else {
        throw std::runtime_error("Duration type not handled");
    }
}

template<typename Scalar, class TimeOrDuration>
constexpr inline Scalar to_sec(const TimeOrDuration& time_or_duration) {
    static_assert(std::is_floating_point_v<TimeOrDuration> || is_time_point_v<TimeOrDuration> ||
                          is_duration_v<TimeOrDuration>,
            "TimeOrDuration must be of floating point scalar, time_point or duration type.");
    if constexpr (std::is_floating_point_v<TimeOrDuration>) {
        return static_cast<Scalar>(time_or_duration);
    } else if constexpr (is_time_point_v<TimeOrDuration>) {
        return to_sec<Scalar, typename TimeOrDuration::duration>(time_or_duration.time_since_epoch());
    } else if constexpr (is_duration_v<TimeOrDuration>) {
        return std::chrono::duration<Scalar>(time_or_duration).count();
    } else {
        throw std::runtime_error("TimeOrDuration type not handled");
    }
}

template<typename Scalar, class TimeOrDuration>
std::vector<Scalar> to_secs(const std::vector<TimeOrDuration>& times_or_durations) {
    std::vector<Scalar> seconds(times_or_durations.size());
    std::transform(times_or_durations.cbegin(), times_or_durations.cend(), seconds.begin(),
            to_sec<Scalar, TimeOrDuration>);
    return seconds;
}

template<class Time, typename Scalar>
constexpr inline Time to_time(const Scalar seconds) {
    static_assert(std::is_floating_point_v<Time> || is_time_point_v<Time>,
            "Time must be of floating point Scalar or time_point type.");
    if constexpr (std::is_floating_point_v<Time>) {
        return static_cast<Time>(seconds);
    } else if constexpr (is_time_point_v<Time>) {
        return Time(to_duration<typename Time::duration>(seconds));
    }
}

template<class Time, typename Scalar = double>
std::vector<Time> to_times(const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& seconds) {
    const std::size_t size = static_cast<std::size_t>(seconds.size());
    std::vector<Time> times(size);
    for (std::size_t i = 0; i < size; ++i) {
        times[i] = to_time<Time, Scalar>(seconds[i]);
    }
    return times;
}

template<class Time, typename Scalar = double>
std::vector<Time> to_times(const Eigen::Matrix<Scalar, 1, Eigen::Dynamic>& seconds) {
    const std::size_t size = static_cast<std::size_t>(seconds.size());
    std::vector<Time> times(size);
    for (std::size_t i = 0; i < size; ++i) {
        times[i] = to_time<Time, Scalar>(seconds[i]);
    }
    return times;
}

}

#endif
