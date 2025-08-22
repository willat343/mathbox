#include "mathbox/time.hpp"

namespace math {

template std::vector<std::chrono::time_point<std::chrono::steady_clock>>
lin_spaced<std::chrono::time_point<std::chrono::steady_clock>, std::chrono::steady_clock::duration>(
        const std::chrono::steady_clock::duration, const std::chrono::time_point<std::chrono::steady_clock>,
        const std::chrono::time_point<std::chrono::steady_clock>);
template std::vector<std::chrono::time_point<std::chrono::system_clock>>
lin_spaced<std::chrono::time_point<std::chrono::system_clock>, std::chrono::system_clock::duration>(
        const std::chrono::system_clock::duration, const std::chrono::time_point<std::chrono::system_clock>,
        const std::chrono::time_point<std::chrono::system_clock>);

template std::vector<std::chrono::time_point<std::chrono::steady_clock>>
range<std::chrono::time_point<std::chrono::steady_clock>, std::chrono::steady_clock::duration>(
        const std::chrono::steady_clock::duration, const std::chrono::time_point<std::chrono::steady_clock>,
        const std::chrono::time_point<std::chrono::steady_clock>);
template std::vector<std::chrono::time_point<std::chrono::system_clock>>
range<std::chrono::time_point<std::chrono::system_clock>, std::chrono::system_clock::duration>(
        const std::chrono::system_clock::duration, const std::chrono::time_point<std::chrono::system_clock>,
        const std::chrono::time_point<std::chrono::system_clock>);

template std::vector<std::chrono::time_point<std::chrono::steady_clock>>
to_times<std::chrono::time_point<std::chrono::steady_clock>, double>(const Eigen::Matrix<double, Eigen::Dynamic, 1>&);
template std::vector<std::chrono::time_point<std::chrono::system_clock>>
to_times<std::chrono::time_point<std::chrono::system_clock>, double>(const Eigen::Matrix<double, Eigen::Dynamic, 1>&);

template std::vector<std::chrono::time_point<std::chrono::steady_clock>>
to_times<std::chrono::time_point<std::chrono::steady_clock>, double>(const Eigen::Matrix<double, 1, Eigen::Dynamic>&);
template std::vector<std::chrono::time_point<std::chrono::system_clock>>
to_times<std::chrono::time_point<std::chrono::system_clock>, double>(const Eigen::Matrix<double, 1, Eigen::Dynamic>&);

}
