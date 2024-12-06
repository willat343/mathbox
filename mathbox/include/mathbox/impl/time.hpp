#ifndef MATHBOX_IMPL_TIME_HPP
#define MATHBOX_IMPL_TIME_HPP

#include <chrono>

#include "mathbox/time.hpp"

namespace math {

template<class Duration, typename Scalar>
Duration to_duration(const Scalar seconds) {
    static_assert(is_duration_v<Duration>, "Duration must be of duration type.");
    return std::chrono::duration_cast<Duration>(std::chrono::duration<Scalar>(seconds));
}

template<class TimeOrDuration, typename Scalar>
Scalar to_sec(const TimeOrDuration& time_or_duration) {
    static_assert(is_time_point_v<TimeOrDuration> || is_duration_v<TimeOrDuration>,
            "TimeOrDuration must be of time_point or duration type.");
    if constexpr (is_time_point_v<TimeOrDuration>) {
        return to_sec(time_or_duration.time_since_epoch());
    } else if constexpr (is_duration_v<TimeOrDuration>) {
        return std::chrono::duration<Scalar>(time_or_duration).count();
    }
}

template<class Time, typename Scalar>
Time to_time(const Scalar seconds) {
    static_assert(is_time_point_v<Time>, "Time must be of time_point type.");
    return Time(to_duration<typename Time::duration>(seconds));
}

}

#endif
