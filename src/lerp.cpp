#include "mathbox/lerp.hpp"

namespace math {

template double lerp<double, double>(const double&, const double&, const double);
template std::chrono::time_point<std::chrono::steady_clock> lerp<std::chrono::time_point<std::chrono::steady_clock>,
        double>(const std::chrono::time_point<std::chrono::steady_clock>&,
        const std::chrono::time_point<std::chrono::steady_clock>&, const double);
template std::chrono::time_point<std::chrono::system_clock> lerp<std::chrono::time_point<std::chrono::system_clock>,
        double>(const std::chrono::time_point<std::chrono::system_clock>&,
        const std::chrono::time_point<std::chrono::system_clock>&, const double);
template std::chrono::nanoseconds lerp<std::chrono::nanoseconds, double>(const std::chrono::nanoseconds&,
        const std::chrono::nanoseconds&, const double);

template double linear_function<double, double>(const double, const double, const double, const double, const double);
template double linear_function<double, std::chrono::time_point<std::chrono::steady_clock>>(
        const std::chrono::time_point<std::chrono::steady_clock>,
        const std::chrono::time_point<std::chrono::steady_clock>,
        const std::chrono::time_point<std::chrono::steady_clock>, const double, const double);
template double linear_function<double, std::chrono::time_point<std::chrono::system_clock>>(
        const std::chrono::time_point<std::chrono::system_clock>,
        const std::chrono::time_point<std::chrono::system_clock>,
        const std::chrono::time_point<std::chrono::system_clock>, const double, const double);

}
