#include "mathbox/vector_operations.hpp"

namespace math {

template Eigen::Matrix<double, Eigen::Dynamic, 1> lin_spaced_vector<double>(const double, const double, const double);

template Eigen::Matrix<double, Eigen::Dynamic, 1> range_vector<double>(const double, const double, const double);

}
