#ifndef MATHBOX_IMPL_NSPHERE_HPP
#define MATHBOX_IMPL_NSPHERE_HPP

#include <Eigen/Core>
#include <cmath>

#include "mathbox/matrix_operations.hpp"
#include "mathbox/nsphere.hpp"

namespace math {

template<int D, typename Scalar>
Eigen::Matrix<Scalar, D, 1> spherical_to_cartesian(const Eigen::Matrix<Scalar, D, 1>& spherical_coordinates) {
    static_assert(D >= 1, "Runtime dimensionality of space in conversion from spherical to cartesian coordinates must "
                          "be at least 1.");
    if constexpr (D == 1) {
        return spherical_coordinates;
    }
    const Eigen::Array<Scalar, D - 1, 1> cos_angles = spherical_coordinates.template tail<D - 1>().array().cos();
    const Eigen::Array<Scalar, D - 1, 1> sin_angles = spherical_coordinates.template tail<D - 1>().array().sin();
    Eigen::Matrix<Scalar, D, 1> cartesian_coordinates;
    cartesian_coordinates[0] = spherical_coordinates[0];
    for (int row = 1; row < D; ++row) {
        cartesian_coordinates[row] = cartesian_coordinates[row - 1] * sin_angles[row - 1];
    }
    for (int row = 0; row < D - 1; ++row) {
        cartesian_coordinates[row] *= cos_angles[row];
    }
    return cartesian_coordinates;
}

template<int D, typename Scalar>
Eigen::Matrix<Scalar, D, 1> cartesian_to_spherical(const Eigen::Matrix<Scalar, D, 1>& cartesian_coordinates) {
    static_assert(D >= 1, "Runtime dimensionality of space in conversion from spherical to cartesian coordinates must "
                          "be at least 1.");
    if constexpr (D == 1) {
        return cartesian_coordinates;
    }
    const Eigen::Matrix<Scalar, D, 1> cartesian_squared = cartesian_coordinates.cwiseProduct(cartesian_coordinates);
    const Eigen::Matrix<Scalar, D, 1> cumulative_cartesian_squared = cumulative_col_top_sum(cartesian_squared);
    const Eigen::Matrix<Scalar, D - 1, 1> sqrt_terms = cumulative_cartesian_squared.template head<D - 1>().cwiseSqrt();
    Eigen::Matrix<Scalar, D, 1> spherical_coordinates;
    spherical_coordinates[0] = sqrt_terms[0];
    for (int row = 1; row < D - 1; ++row) {
        // Since the sqrt terms are >= 0, the result of atan2 will be in desired range [0, pi].
        spherical_coordinates[row] = std::atan2(sqrt_terms[row], cartesian_coordinates[row - 1]);
    }
    // Last term does not use the square root term. We leave it in range [-pi, pi].
    spherical_coordinates[D - 1] = std::atan2(cartesian_coordinates[D - 1], cartesian_coordinates[D - 2]);
    return spherical_coordinates;
}

template<int D, typename Scalar>
Eigen::Matrix<Scalar, D - 1, 1> cartesian_to_spherical_angles(
        const Eigen::Matrix<Scalar, D, 1>& unit_cartesian_coordinates) {
    static_assert(D >= 2, "Runtime dimensionality of space in conversion from unit cartesian coordinates to spherical "
                          "angles must be at least 2.");
    const Eigen::Matrix<Scalar, D - 1, 1> cartesian_squared =
            unit_cartesian_coordinates.template tail<D - 1>().cwiseProduct(
                    unit_cartesian_coordinates.template tail<D - 1>());
    const Eigen::Matrix<Scalar, D - 1, 1> cumulative_cartesian_squared = cumulative_col_top_sum(cartesian_squared);
    const Eigen::Matrix<Scalar, D - 2, 1> sqrt_terms = cumulative_cartesian_squared.template head<D - 2>().cwiseSqrt();
    Eigen::Matrix<Scalar, D - 1, 1> spherical_angles;
    for (int row = 0; row < D - 2; ++row) {
        // Since the sqrt terms are >= 0, the result of atan2 will be in desired range [0, pi].
        spherical_angles[row] = std::atan2(sqrt_terms[row], unit_cartesian_coordinates[row]);
    }
    // Last term does not use the square root term. We leave it in range [-pi, pi].
    spherical_angles[D - 2] = std::atan2(unit_cartesian_coordinates[D - 1], unit_cartesian_coordinates[D - 2]);
    return spherical_angles;
}

template<int D, typename Scalar>
inline Eigen::Matrix<Scalar, D, 1> spherical_angles_to_unit_cartesian(
        const Eigen::Matrix<Scalar, D - 1, 1>& spherical_angles) {
    static_assert(D >= 2, "Runtime dimensionality of space in conversion from spherical angles to unit cartesian "
                          "coordinates must be at least 2.");
    Eigen::Matrix<Scalar, D, 1> spherical_coordinates;
    spherical_coordinates << static_cast<Scalar>(1.0), spherical_angles;
    return spherical_to_cartesian<D, Scalar>(spherical_coordinates);
}

}

#if !MATHBOX_HEADER_ONLY
namespace math {

extern template Eigen::Matrix<double, 1, 1> cartesian_to_spherical<1, double>(const Eigen::Matrix<double, 1, 1>&);
extern template Eigen::Matrix<double, 2, 1> cartesian_to_spherical<2, double>(const Eigen::Matrix<double, 2, 1>&);
extern template Eigen::Matrix<double, 3, 1> cartesian_to_spherical<3, double>(const Eigen::Matrix<double, 3, 1>&);

extern template Eigen::Matrix<double, 1, 1> cartesian_to_spherical_angles<2, double>(
        const Eigen::Matrix<double, 2, 1>&);
extern template Eigen::Matrix<double, 2, 1> cartesian_to_spherical_angles<3, double>(
        const Eigen::Matrix<double, 3, 1>&);

extern template Eigen::Matrix<double, 2, 1> spherical_angles_to_unit_cartesian<2, double>(
        const Eigen::Matrix<double, 1, 1>&);
extern template Eigen::Matrix<double, 3, 1> spherical_angles_to_unit_cartesian<3, double>(
        const Eigen::Matrix<double, 2, 1>&);

extern template Eigen::Matrix<double, 1, 1> spherical_to_cartesian<1, double>(const Eigen::Matrix<double, 1, 1>&);
extern template Eigen::Matrix<double, 2, 1> spherical_to_cartesian<2, double>(const Eigen::Matrix<double, 2, 1>&);
extern template Eigen::Matrix<double, 3, 1> spherical_to_cartesian<3, double>(const Eigen::Matrix<double, 3, 1>&);

}
#endif

#endif
