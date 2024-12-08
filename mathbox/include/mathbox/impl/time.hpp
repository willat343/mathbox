#ifndef MATHBOX_IMPL_TIME_HPP
#define MATHBOX_IMPL_TIME_HPP

#include <chrono>

#include "mathbox/time.hpp"

namespace math {

template<typename Scalar, class Duration>
Scalar fraction(const Duration numerator, const Duration denominator) {
    return to_sec<Scalar>(numerator) / to_sec<Scalar>(denominator);
}

template<class Duration, typename Scalar>
Duration to_duration(const Scalar seconds) {
    static_assert(std::is_floating_point_v<Duration> || is_duration_v<Duration>,
            "Duration must be of floating point scalar or duration type.");
    if constexpr (std::is_floating_point_v<Duration>) {
        return static_cast<Duration>(seconds);
    } else if constexpr (is_duration_v<Duration>) {
        return std::chrono::duration_cast<Duration>(std::chrono::duration<Scalar>(seconds));
    }
}

template<typename Scalar, class TimeOrDuration>
Scalar to_sec(const TimeOrDuration& time_or_duration) {
    static_assert(std::is_floating_point_v<TimeOrDuration> || is_time_point_v<TimeOrDuration> ||
                          is_duration_v<TimeOrDuration>,
            "TimeOrDuration must be of floating point scalar, time_point or duration type.");
    if constexpr (std::is_floating_point_v<TimeOrDuration>) {
        return static_cast<Scalar>(time_or_duration);
    } else if constexpr (is_time_point_v<TimeOrDuration>) {
        return to_sec<Scalar, typename TimeOrDuration::duration>(time_or_duration.time_since_epoch());
    } else if constexpr (is_duration_v<TimeOrDuration>) {
        return std::chrono::duration<Scalar>(time_or_duration).count();
    }
}

template<class Time, typename Scalar>
Time to_time(const Scalar seconds) {
    static_assert(std::is_floating_point_v<Time> || is_time_point_v<Time>,
            "Time must be of floating point Scalar or time_point type.");
    if constexpr (std::is_floating_point_v<Time>) {
        return static_cast<Time>(seconds);
    } else if constexpr (is_time_point_v<Time>) {
        return Time(to_duration<typename Time::duration>(seconds));
    }
}

}

#endif
