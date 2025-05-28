#ifndef MATHBOX_IMPL_TIME_HPP
#define MATHBOX_IMPL_TIME_HPP

#include "mathbox/time.hpp"
#include "mathbox/vector_operations.hpp"

namespace math {

template<class Time, class Duration>
std::vector<Time> lin_spaced(const Duration step, const Time start, const Time end) {
    Eigen::Matrix<double, Eigen::Dynamic, 1> lin_spaced_vector =
            lin_spaced<double>(cppbox::to_sec(step), cppbox::to_sec(start), cppbox::to_sec(end));
    const std::size_t size = static_cast<std::size_t>(lin_spaced_vector.size());
    std::vector<Time> lin_spaced_times(size);
    // Eigen 3.3 does not support iterators so cannot use std::transform.
    for (std::size_t i = 0; i < size; ++i) {
        lin_spaced_times[i] = cppbox::to_time<Time>(lin_spaced_vector[i]);
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

template<class Time, typename Scalar = double>
std::vector<Time> to_times(const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& seconds) {
    const std::size_t size = static_cast<std::size_t>(seconds.size());
    std::vector<Time> times(size);
    for (std::size_t i = 0; i < size; ++i) {
        times[i] = cppbox::to_time<Time, Scalar>(seconds[i]);
    }
    return times;
}

template<class Time, typename Scalar = double>
std::vector<Time> to_times(const Eigen::Matrix<Scalar, 1, Eigen::Dynamic>& seconds) {
    const std::size_t size = static_cast<std::size_t>(seconds.size());
    std::vector<Time> times(size);
    for (std::size_t i = 0; i < size; ++i) {
        times[i] = cppbox::to_time<Time, Scalar>(seconds[i]);
    }
    return times;
}

}

#endif
