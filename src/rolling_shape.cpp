#include "mathbox/rolling_shape.hpp"

namespace math {

template class RollingShape<double>;
template class Sphere<double>;
template class Torus<double>;

template Eigen::Vector3d Sphere<double>::contact_to_centre<double>(const Eigen::Quaterniond&,
        const Eigen::Vector3d&) const;

template Eigen::Vector3d Torus<double>::contact_to_centre<double>(const Eigen::Quaterniond&,
        const Eigen::Vector3d&) const;

template Eigen::Vector3d rolling_displacement<double, Sphere<double>>(const Sphere<double>&, const Eigen::Quaterniond&,
        const Eigen::Quaterniond&, const Eigen::Vector3d&);
template Eigen::Vector3d rolling_displacement<double, Torus<double>>(const Torus<double>&, const Eigen::Quaterniond&,
        const Eigen::Quaterniond&, const Eigen::Vector3d&);

template Eigen::Vector3d rolling_displacement<double, Sphere<double>, std::vector<Eigen::Quaterniond>>(
        const Sphere<double>&, const std::vector<Eigen::Quaterniond>&, const Eigen::Vector3d&);
template Eigen::Vector3d rolling_displacement<double, Torus<double>, std::vector<Eigen::Quaterniond>>(
        const Torus<double>&, const std::vector<Eigen::Quaterniond>&, const Eigen::Vector3d&);

}
