#include "mathbox/nsphere.hpp"

namespace math {

template Eigen::Matrix<double, 1, 1> cartesian_to_spherical<1, double>(const Eigen::Matrix<double, 1, 1>&);
template Eigen::Matrix<double, 2, 1> cartesian_to_spherical<2, double>(const Eigen::Matrix<double, 2, 1>&);
template Eigen::Matrix<double, 3, 1> cartesian_to_spherical<3, double>(const Eigen::Matrix<double, 3, 1>&);

template Eigen::Matrix<double, 1, 1> cartesian_to_spherical_angles<2, double>(const Eigen::Matrix<double, 2, 1>&);
template Eigen::Matrix<double, 2, 1> cartesian_to_spherical_angles<3, double>(const Eigen::Matrix<double, 3, 1>&);

template Eigen::Matrix<double, 2, 1> spherical_angles_to_unit_cartesian<2, double>(const Eigen::Matrix<double, 1, 1>&);
template Eigen::Matrix<double, 3, 1> spherical_angles_to_unit_cartesian<3, double>(const Eigen::Matrix<double, 2, 1>&);

template Eigen::Matrix<double, 1, 1> spherical_to_cartesian<1, double>(const Eigen::Matrix<double, 1, 1>&);
template Eigen::Matrix<double, 2, 1> spherical_to_cartesian<2, double>(const Eigen::Matrix<double, 2, 1>&);
template Eigen::Matrix<double, 3, 1> spherical_to_cartesian<3, double>(const Eigen::Matrix<double, 3, 1>&);

}
