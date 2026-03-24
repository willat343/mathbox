#ifndef MATHBOX_IMPL_TIME_HPP
#define MATHBOX_IMPL_TIME_HPP

#include <cppbox/exceptions.hpp>

#include "mathbox/time.hpp"
#include "mathbox/vector_operations.hpp"

namespace math {

template<class Time, class Duration>
std::vector<Time> lin_spaced(const Duration step, const Time start, const Time end) {
    Eigen::Matrix<double, Eigen::Dynamic, 1> vector =
            lin_spaced_vector<double>(cppbox::to_sec(step), cppbox::to_sec(start), cppbox::to_sec(end));
    const std::size_t size = static_cast<std::size_t>(vector.size());
    std::vector<Time> times(size);
    // Eigen 3.3 does not support iterators so cannot use std::transform.
    for (std::size_t i = 0; i < size; ++i) {
        times[i] = cppbox::to_time<Time>(vector[i]);
    }
    return times;
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

template<class Time, typename Scalar = double>
Eigen::Matrix<Scalar, Eigen::Dynamic, 1> from_times(const std::vector<Time>& times) {
    const int size = static_cast<int>(times.size());
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> vector(size, 1);
    for (int i = 0; i < size; ++i) {
        vector[i] = cppbox::to_sec(times[i]);
    }
    return vector;
}

template<class Time, typename Scalar = double>
Eigen::Matrix<Scalar, Eigen::Dynamic, 1> from_times(const std::vector<Time>& times, const int start_index) {
    throw_if(start_index > static_cast<int>(times.size()), "Requested times outside of bounds.");
    const int size = static_cast<int>(times.size()) - start_index;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> vector(size, 1);
    for (int i = 0; i < size; ++i) {
        vector[i] = cppbox::to_sec(times[start_index + i]);
    }
    return vector;
}

template<class Time, typename Scalar = double>
Eigen::Matrix<Scalar, Eigen::Dynamic, 1> from_times(const std::vector<Time>& times, const int start_index,
        const int size) {
    throw_if(start_index + size > static_cast<int>(times.size()), "Requested times outside of bounds.");
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> vector(size, 1);
    for (int i = 0; i < size; ++i) {
        vector[i] = cppbox::to_sec(times[start_index + i]);
    }
    return vector;
}

template<class Time, int Size, typename Scalar = double>
Eigen::Matrix<Scalar, Size, 1> from_times(const std::vector<Time>& times) {
    Eigen::Matrix<Scalar, Size, 1> vector;
    for (int i = 0; i < Size; ++i) {
        vector[i] = cppbox::to_sec(times[i]);
    }
    return vector;
}

template<class Time, int Size, typename Scalar = double>
Eigen::Matrix<Scalar, Size, 1> from_times(const std::vector<Time>& times, const int start_index) {
    throw_if(start_index + Size > static_cast<int>(times.size()), "Requested times outside of bounds.");
    Eigen::Matrix<Scalar, Size, 1> vector;
    for (int i = 0; i < Size; ++i) {
        vector[i] = cppbox::to_sec(times[start_index + i]);
    }
    return vector;
}

}

#if !MATHBOX_HEADER_ONLY
namespace math {

#ifndef DOXYGEN_EXCLUDE
extern template std::vector<std::chrono::time_point<std::chrono::steady_clock>>
lin_spaced<std::chrono::time_point<std::chrono::steady_clock>, std::chrono::steady_clock::duration>(
        const std::chrono::steady_clock::duration, const std::chrono::time_point<std::chrono::steady_clock>,
        const std::chrono::time_point<std::chrono::steady_clock>);
extern template std::vector<std::chrono::time_point<std::chrono::system_clock>>
lin_spaced<std::chrono::time_point<std::chrono::system_clock>, std::chrono::system_clock::duration>(
        const std::chrono::system_clock::duration, const std::chrono::time_point<std::chrono::system_clock>,
        const std::chrono::time_point<std::chrono::system_clock>);
#endif

#ifndef DOXYGEN_EXCLUDE
extern template std::vector<std::chrono::time_point<std::chrono::steady_clock>>
range<std::chrono::time_point<std::chrono::steady_clock>, std::chrono::steady_clock::duration>(
        const std::chrono::steady_clock::duration, const std::chrono::time_point<std::chrono::steady_clock>,
        const std::chrono::time_point<std::chrono::steady_clock>);
extern template std::vector<std::chrono::time_point<std::chrono::system_clock>>
range<std::chrono::time_point<std::chrono::system_clock>, std::chrono::system_clock::duration>(
        const std::chrono::system_clock::duration, const std::chrono::time_point<std::chrono::system_clock>,
        const std::chrono::time_point<std::chrono::system_clock>);
#endif

extern template std::vector<std::chrono::time_point<std::chrono::steady_clock>>
to_times<std::chrono::time_point<std::chrono::steady_clock>, double>(const Eigen::Matrix<double, Eigen::Dynamic, 1>&);
extern template std::vector<std::chrono::time_point<std::chrono::system_clock>>
to_times<std::chrono::time_point<std::chrono::system_clock>, double>(const Eigen::Matrix<double, Eigen::Dynamic, 1>&);

extern template std::vector<std::chrono::time_point<std::chrono::steady_clock>>
to_times<std::chrono::time_point<std::chrono::steady_clock>, double>(const Eigen::Matrix<double, 1, Eigen::Dynamic>&);
extern template std::vector<std::chrono::time_point<std::chrono::system_clock>>
to_times<std::chrono::time_point<std::chrono::system_clock>, double>(const Eigen::Matrix<double, 1, Eigen::Dynamic>&);

}
#endif

#endif
