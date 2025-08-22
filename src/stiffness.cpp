#include "mathbox/stiffness.hpp"

namespace math {

template Eigen::MatrixXd stiffness_from_sigmas<Eigen::VectorXd>(const Eigen::MatrixBase<Eigen::VectorXd>&);

template Eigen::MatrixXd stiffness_from_sigmas<Eigen::Dynamic>(const Eigen::Ref<const Eigen::VectorXd>&);

template double stiffness_from_sigma<double>(const double);

template Eigen::MatrixXd stiffness_from_variances<Eigen::VectorXd>(const Eigen::MatrixBase<Eigen::VectorXd>&);

template Eigen::MatrixXd stiffness_from_variances<Eigen::Dynamic>(const Eigen::Ref<const Eigen::VectorXd>&);

template double stiffness_from_variance<double>(const double);

}
