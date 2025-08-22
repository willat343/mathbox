#include "mathbox/covariance.hpp"

namespace math {

template class PositiveSemiDefiniteMatrix<double, Eigen::Dynamic>;
template class PositiveSemiDefiniteMatrix<double, 1>;
template class PositiveSemiDefiniteMatrix<double, 2>;
template class PositiveSemiDefiniteMatrix<double, 3>;
template class PositiveSemiDefiniteMatrix<double, 4>;
template class PositiveSemiDefiniteMatrix<double, 5>;
template class PositiveSemiDefiniteMatrix<double, 6>;
template class PositiveSemiDefiniteMatrix<double, 7>;
template class PositiveSemiDefiniteMatrix<double, 8>;
template class PositiveSemiDefiniteMatrix<double, 9>;

template class CovarianceDensity<double, Eigen::Dynamic>;
template class CovarianceDensity<double, 1>;
template class CovarianceDensity<double, 2>;
template class CovarianceDensity<double, 3>;
template class CovarianceDensity<double, 4>;
template class CovarianceDensity<double, 5>;
template class CovarianceDensity<double, 6>;
template class CovarianceDensity<double, 7>;
template class CovarianceDensity<double, 8>;
template class CovarianceDensity<double, 9>;

}
