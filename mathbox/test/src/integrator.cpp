#include "mathbox/integrator.hpp"

#include <gtest/gtest.h>

#include <chrono>

template<typename MathType>
MathType zero_instance() {
    static_assert(math::is_math_type_v<MathType>, "MathType was not a math type.");
    if constexpr (std::is_base_of_v<Eigen::MatrixBase<MathType>, MathType>) {
        // Set any dynamic-sized dimensions to 1
        return math::MathTypeTraits<MathType>::zero(
                MathType::RowsAtCompileTime == Eigen::Dynamic ? 1 : MathType::RowsAtCompileTime,
                MathType::ColsAtCompileTime == Eigen::Dynamic ? 1 : MathType::ColsAtCompileTime);
    } else {
        return math::MathTypeTraits<MathType>::zero();
    }
}

template<typename Integrator>
void test_integrate_constant(const typename Integrator::IndependentVariableType& start,
        const typename Integrator::IndependentVariableType& end, const typename Integrator::MathType constant,
        const typename Integrator::IndependentVariableDifferenceType integration_step,
        const typename Integrator::MathType& initial_value,
        const typename Integrator::ArithmeticTypeScalar precision_per_unit_area) {
    using MathType = Integrator::MathType;
    using IndependentVariableType = Integrator::IndependentVariableType;
    using ArithmeticTypeScalar = Integrator::ArithmeticTypeScalar;
    Integrator integrator;
    MathType expected_value = constant;
    const ArithmeticTypeScalar duration = math::to_sec(end - start);
    expected_value *= duration;
    expected_value += initial_value;
    const MathType value = integrator.integrate(start, end, integration_step, initial_value,
            [constant](const IndependentVariableType) -> MathType { return constant; });
    if constexpr (std::is_floating_point_v<MathType>) {
        const ArithmeticTypeScalar precision = precision_per_unit_area * std::abs(duration) * std::abs(constant);
        EXPECT_NEAR(value, expected_value, precision);
    } else if constexpr (std::is_base_of_v<Eigen::MatrixBase<MathType>, MathType>) {
        const ArithmeticTypeScalar precision =
                precision_per_unit_area * std::abs(duration) * std::abs(constant.value());
        // Note isApprox does not work when integration is close to zero
        EXPECT_NEAR(value.value(), expected_value.value(), precision);
    }
}

template<typename Integrator>
void test_integrate_linear(const typename Integrator::IndependentVariableType& start,
        const typename Integrator::IndependentVariableType& end, const typename Integrator::MathType start_height,
        const typename Integrator::MathType end_height,
        const typename Integrator::IndependentVariableDifferenceType integration_step,
        const typename Integrator::MathType& initial_value,
        const typename Integrator::ArithmeticTypeScalar precision_per_unit_area) {
    using MathType = Integrator::MathType;
    using IndependentVariableType = Integrator::IndependentVariableType;
    using ArithmeticTypeScalar = Integrator::ArithmeticTypeScalar;
    Integrator integrator;
    const MathType mean_height = (start_height + end_height) / ArithmeticTypeScalar(2.0);
    MathType expected_value = mean_height;
    const ArithmeticTypeScalar duration = math::to_sec(end - start);
    expected_value *= duration;
    expected_value += initial_value;
    const MathType value1 = integrator.integrate(start, end, integration_step, initial_value,
            std::bind(math::linear_function<MathType, IndependentVariableType>, std::placeholders::_1, start, end,
                    start_height, end_height));
    const MathType value2 = integrator.integrate(start, end, integration_step, initial_value,
            [start, end, start_height, end_height](const IndependentVariableType t) -> MathType {
                return math::linear_function(t, start, end, start_height, end_height);
            });
    if constexpr (std::is_floating_point_v<MathType>) {
        // Absolute area is more complicated than integral, depending on if a zero-crossing occurred
        ArithmeticTypeScalar area;
        if (start_height * end_height >= ArithmeticTypeScalar(0.0)) {
            area = std::abs(duration) * std::abs(mean_height);
        } else {
            const ArithmeticTypeScalar fraction_before_crossing = std::abs(start_height / (start_height - end_height));
            area = std::abs(duration) *
                   (fraction_before_crossing * std::abs(start_height) +
                           (ArithmeticTypeScalar(1.0) - fraction_before_crossing) * std::abs(end_height)) /
                   ArithmeticTypeScalar(2.0);
        }
        const ArithmeticTypeScalar precision = precision_per_unit_area * area;
        EXPECT_EQ(value1, value2);
        EXPECT_NEAR(value1, expected_value, precision);
        EXPECT_NEAR(value2, expected_value, precision);
    } else if constexpr (std::is_base_of_v<Eigen::MatrixBase<MathType>, MathType>) {
        // Absolute area is more complicated than integral, depending on if a zero-crossing occurred
        ArithmeticTypeScalar area;
        if (start_height.value() * end_height.value() >= ArithmeticTypeScalar(0.0)) {
            area = std::abs(duration) * std::abs(mean_height.value());
        } else {
            const ArithmeticTypeScalar fraction_before_crossing =
                    std::abs(start_height.value() / (start_height.value() - end_height.value()));
            area = std::abs(duration) *
                   (fraction_before_crossing * std::abs(start_height.value()) +
                           (ArithmeticTypeScalar(1.0) - fraction_before_crossing) * std::abs(end_height.value())) /
                   ArithmeticTypeScalar(2.0);
        }
        const ArithmeticTypeScalar precision = precision_per_unit_area * area;
        // Note isApprox does not work when integration is close to zero
        EXPECT_EQ(value1.value(), value2.value());
        EXPECT_NEAR(value1.value(), expected_value.value(), precision);
        EXPECT_NEAR(value2.value(), expected_value.value(), precision);
    }
}

#define GENERATE_INTEGRATION_TEST_CONSTANT_DOUBLE(name, test_postfix, start, end, constant, integration_step,          \
        precision_per_unit_area)                                                                                       \
                                                                                                                       \
    TEST(name##_integration, integrate_constant_double_double_##test_postfix) {                                        \
        test_integrate_constant<math::newton_cotes::name::Integrator<double, double>>(start, end, constant,            \
                math::to_duration<double>(integration_step), zero_instance<double>(), precision_per_unit_area);        \
    }                                                                                                                  \
                                                                                                                       \
    TEST(name##_integration, integrate_constant_fixedmatrix_double_##test_postfix) {                                   \
        test_integrate_constant<math::newton_cotes::name::Integrator<Eigen::Matrix<double, 1, 1>, double>>(start, end, \
                Eigen::Matrix<double, 1, 1>::Constant(constant), math::to_duration<double>(integration_step),          \
                zero_instance<Eigen::Matrix<double, 1, 1>>(), precision_per_unit_area);                                \
    }                                                                                                                  \
                                                                                                                       \
    TEST(name##_integration, integrate_constant_dynamicvector_double_##test_postfix) {                                 \
        test_integrate_constant<                                                                                       \
                math::newton_cotes::name::Integrator<Eigen::Matrix<double, Eigen::Dynamic, 1>, double>>(start, end,    \
                Eigen::Matrix<double, Eigen::Dynamic, 1>::Constant(1, constant),                                       \
                math::to_duration<double>(integration_step),                                                           \
                zero_instance<Eigen::Matrix<double, Eigen::Dynamic, 1>>(), precision_per_unit_area);                   \
    }                                                                                                                  \
                                                                                                                       \
    TEST(name##_integration, integrate_constant_dynamicmatrix_double_##test_postfix) {                                 \
        test_integrate_constant<                                                                                       \
                math::newton_cotes::name::Integrator<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>, double>>(  \
                start, end, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Constant(1, 1, constant),           \
                math::to_duration<double>(integration_step),                                                           \
                zero_instance<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(), precision_per_unit_area);      \
    }

#define GENERATE_INTEGRATION_TEST_CONSTANT_TIME(name, test_postfix, start, end, constant, integration_step,        \
        precision_per_unit_area)                                                                                   \
                                                                                                                   \
    TEST(name##_integration, integrate_constant_double_time_##test_postfix) {                                      \
        test_integrate_constant<                                                                                   \
                math::newton_cotes::name::Integrator<double, std::chrono::time_point<std::chrono::steady_clock>>>( \
                math::to_time<std::chrono::time_point<std::chrono::steady_clock>>(start),                          \
                math::to_time<std::chrono::time_point<std::chrono::steady_clock>>(end), constant,                  \
                math::to_duration<std::chrono::time_point<std::chrono::steady_clock>::duration>(integration_step), \
                zero_instance<double>(), precision_per_unit_area);                                                 \
    }                                                                                                              \
                                                                                                                   \
    TEST(name##_integration, integrate_constant_fixedmatrix_time_##test_postfix) {                                 \
        test_integrate_constant<math::newton_cotes::name::Integrator<Eigen::Matrix<double, 1, 1>,                  \
                std::chrono::time_point<std::chrono::steady_clock>>>(                                              \
                math::to_time<std::chrono::time_point<std::chrono::steady_clock>>(start),                          \
                math::to_time<std::chrono::time_point<std::chrono::steady_clock>>(end),                            \
                Eigen::Matrix<double, 1, 1>::Constant(constant),                                                   \
                math::to_duration<std::chrono::time_point<std::chrono::steady_clock>::duration>(integration_step), \
                zero_instance<Eigen::Matrix<double, 1, 1>>(), precision_per_unit_area);                            \
    }                                                                                                              \
                                                                                                                   \
    TEST(name##_integration, integrate_constant_dynamicvector_time_##test_postfix) {                               \
        test_integrate_constant<math::newton_cotes::name::Integrator<Eigen::Matrix<double, Eigen::Dynamic, 1>,     \
                std::chrono::time_point<std::chrono::steady_clock>>>(                                              \
                math::to_time<std::chrono::time_point<std::chrono::steady_clock>>(start),                          \
                math::to_time<std::chrono::time_point<std::chrono::steady_clock>>(end),                            \
                Eigen::Matrix<double, Eigen::Dynamic, 1>::Constant(1, constant),                                   \
                math::to_duration<std::chrono::time_point<std::chrono::steady_clock>::duration>(integration_step), \
                zero_instance<Eigen::Matrix<double, Eigen::Dynamic, 1>>(), precision_per_unit_area);               \
    }                                                                                                              \
                                                                                                                   \
    TEST(name##_integration, integrate_constant_dynamicmatrix_time_##test_postfix) {                               \
        test_integrate_constant<                                                                                   \
                math::newton_cotes::name::Integrator<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>,        \
                        std::chrono::time_point<std::chrono::steady_clock>>>(                                      \
                math::to_time<std::chrono::time_point<std::chrono::steady_clock>>(start),                          \
                math::to_time<std::chrono::time_point<std::chrono::steady_clock>>(end),                            \
                Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Constant(1, 1, constant),                   \
                math::to_duration<std::chrono::time_point<std::chrono::steady_clock>::duration>(integration_step), \
                zero_instance<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(), precision_per_unit_area);  \
    }

#define GENERATE_INTEGRATION_TEST_LINEAR_DOUBLE(name, test_postfix, start, end, start_height, end_height,              \
        integration_step, precision_per_unit_area)                                                                     \
                                                                                                                       \
    TEST(name##_integration, integrate_linear_double_double_##test_postfix) {                                          \
        test_integrate_linear<math::newton_cotes::name::Integrator<double, double>>(start, end, start_height,          \
                end_height, math::to_duration<double>(integration_step), zero_instance<double>(),                      \
                precision_per_unit_area);                                                                              \
    }                                                                                                                  \
                                                                                                                       \
    TEST(name##_integration, integrate_linear_fixedmatrix_double_##test_postfix) {                                     \
        test_integrate_linear<math::newton_cotes::name::Integrator<Eigen::Matrix<double, 1, 1>, double>>(start, end,   \
                Eigen::Matrix<double, 1, 1>::Constant(start_height),                                                   \
                Eigen::Matrix<double, 1, 1>::Constant(end_height), math::to_duration<double>(integration_step),        \
                zero_instance<Eigen::Matrix<double, 1, 1>>(), precision_per_unit_area);                                \
    }                                                                                                                  \
                                                                                                                       \
    TEST(name##_integration, integrate_linear_dynamicvector_double_##test_postfix) {                                   \
        test_integrate_linear<math::newton_cotes::name::Integrator<Eigen::Matrix<double, Eigen::Dynamic, 1>, double>>( \
                start, end, Eigen::Matrix<double, Eigen::Dynamic, 1>::Constant(1, start_height),                       \
                Eigen::Matrix<double, Eigen::Dynamic, 1>::Constant(1, end_height),                                     \
                math::to_duration<double>(integration_step),                                                           \
                zero_instance<Eigen::Matrix<double, Eigen::Dynamic, 1>>(), precision_per_unit_area);                   \
    }                                                                                                                  \
                                                                                                                       \
    TEST(name##_integration, integrate_linear_dynamicmatrix_double_##test_postfix) {                                   \
        test_integrate_linear<                                                                                         \
                math::newton_cotes::name::Integrator<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>, double>>(  \
                start, end, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Constant(1, 1, start_height),       \
                Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Constant(1, 1, end_height),                     \
                math::to_duration<double>(integration_step),                                                           \
                zero_instance<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(), precision_per_unit_area);      \
    }

#define GENERATE_INTEGRATION_TEST_LINEAR_TIME(name, test_postfix, start, end, start_height, end_height,            \
        integration_step, precision_per_unit_area)                                                                 \
                                                                                                                   \
    TEST(name##_integration, integrate_linear_double_time_##test_postfix) {                                        \
        test_integrate_linear<                                                                                     \
                math::newton_cotes::name::Integrator<double, std::chrono::time_point<std::chrono::steady_clock>>>( \
                math::to_time<std::chrono::time_point<std::chrono::steady_clock>>(start),                          \
                math::to_time<std::chrono::time_point<std::chrono::steady_clock>>(end), start_height, end_height,  \
                math::to_duration<std::chrono::time_point<std::chrono::steady_clock>::duration>(integration_step), \
                zero_instance<double>(), precision_per_unit_area);                                                 \
    }                                                                                                              \
                                                                                                                   \
    TEST(name##_integration, integrate_linear_fixedmatrix_time_##test_postfix) {                                   \
        test_integrate_linear<math::newton_cotes::name::Integrator<Eigen::Matrix<double, 1, 1>,                    \
                std::chrono::time_point<std::chrono::steady_clock>>>(                                              \
                math::to_time<std::chrono::time_point<std::chrono::steady_clock>>(start),                          \
                math::to_time<std::chrono::time_point<std::chrono::steady_clock>>(end),                            \
                Eigen::Matrix<double, 1, 1>::Constant(start_height),                                               \
                Eigen::Matrix<double, 1, 1>::Constant(end_height),                                                 \
                math::to_duration<std::chrono::time_point<std::chrono::steady_clock>::duration>(integration_step), \
                zero_instance<Eigen::Matrix<double, 1, 1>>(), precision_per_unit_area);                            \
    }                                                                                                              \
                                                                                                                   \
    TEST(name##_integration, integrate_linear_dynamicvector_time_##test_postfix) {                                 \
        test_integrate_linear<math::newton_cotes::name::Integrator<Eigen::Matrix<double, Eigen::Dynamic, 1>,       \
                std::chrono::time_point<std::chrono::steady_clock>>>(                                              \
                math::to_time<std::chrono::time_point<std::chrono::steady_clock>>(start),                          \
                math::to_time<std::chrono::time_point<std::chrono::steady_clock>>(end),                            \
                Eigen::Matrix<double, Eigen::Dynamic, 1>::Constant(1, start_height),                               \
                Eigen::Matrix<double, Eigen::Dynamic, 1>::Constant(1, end_height),                                 \
                math::to_duration<std::chrono::time_point<std::chrono::steady_clock>::duration>(integration_step), \
                zero_instance<Eigen::Matrix<double, Eigen::Dynamic, 1>>(), precision_per_unit_area);               \
    }                                                                                                              \
                                                                                                                   \
    TEST(name##_integration, integrate_linear_dynamicmatrix_time_##test_postfix) {                                 \
        test_integrate_linear<                                                                                     \
                math::newton_cotes::name::Integrator<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>,        \
                        std::chrono::time_point<std::chrono::steady_clock>>>(                                      \
                math::to_time<std::chrono::time_point<std::chrono::steady_clock>>(start),                          \
                math::to_time<std::chrono::time_point<std::chrono::steady_clock>>(end),                            \
                Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Constant(1, 1, start_height),               \
                Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Constant(1, 1, end_height),                 \
                math::to_duration<std::chrono::time_point<std::chrono::steady_clock>::duration>(integration_step), \
                zero_instance<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(), precision_per_unit_area);  \
    }

#define GENERATE_INTEGRATION_TEST_CONSTANT_STEP(name, step_name, test_postfix, start, end, constant, step,  \
        precision_per_unit_area, time_precision_factor)                                                     \
    GENERATE_INTEGRATION_TEST_CONSTANT_DOUBLE(name, step_name##_##test_postfix, start, end, constant, step, \
            precision_per_unit_area)                                                                        \
    GENERATE_INTEGRATION_TEST_CONSTANT_TIME(name, step_name##_##test_postfix, start, end, constant, step,   \
            time_precision_factor* precision_per_unit_area)

#define GENERATE_INTEGRATION_TEST_LINEAR_STEP(name, step_name, test_postfix, start, end, start_height, end_height,  \
        step, precision_per_unit_area, time_precision_factor)                                                       \
    GENERATE_INTEGRATION_TEST_LINEAR_DOUBLE(name, step_name##_##test_postfix, start, end, start_height, end_height, \
            step, precision_per_unit_area)                                                                          \
    GENERATE_INTEGRATION_TEST_LINEAR_TIME(name, step_name##_##test_postfix, start, end, start_height, end_height,   \
            step, time_precision_factor* precision_per_unit_area)

GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, us, v0, 0.0, 1.0, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, us, v1, 0.0, 10.0, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, us, v2, 0.0, 1.0, 10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, us, v3, 0.0, 1.0, 100.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, us, v4, -1.0, 1.0, 10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, us, v5, 0.123456789, 0.987654321, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, us, v6, 1.0, -1.0, 10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, us, v7, 0.0, 1.0, -1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, us, v8, 0.0, 1.0, -10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, us, v9, 0.0, 1.0, -100.0, 1.0e-6, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, us, v0, 0.0, 1.0, 1.0, 2.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, us, v1, 0.0, 10.0, 1.0, 2.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, us, v2, 0.0, 1.0, 10.0, 20.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, us, v3, 0.0, 1.0, 100.0, 200.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, us, v4, -1.0, 1.0, 10.0, 0.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, us, v5, 0.123456789, 0.987654321, 1.0, 2.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, us, v6, 1.0, -1.0, 10.0, 20.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, us, v7, 0.0, 1.0, -1.0, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, us, v8, 0.0, 1.0, 10.0, -10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, us, v9, 0.0, 1.0, -100.0, 100.0, 1.0e-6, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, ms, v0, 0.0, 1.0, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, ms, v1, 0.0, 10.0, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, ms, v2, 0.0, 1.0, 10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, ms, v3, 0.0, 1.0, 100.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, ms, v4, -1.0, 1.0, 10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, ms, v5, 0.123456789, 0.987654321, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, ms, v6, 1.0, -1.0, 10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, ms, v7, 0.0, 1.0, -1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, ms, v8, 0.0, 1.0, -10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, ms, v9, 0.0, 1.0, -100.0, 1.0e-3, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, ms, v0, 0.0, 1.0, 1.0, 2.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, ms, v1, 0.0, 10.0, 1.0, 2.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, ms, v2, 0.0, 1.0, 10.0, 20.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, ms, v3, 0.0, 1.0, 100.0, 200.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, ms, v4, -1.0, 1.0, 10.0, 0.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, ms, v5, 0.123456789, 0.987654321, 1.0, 2.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, ms, v6, 1.0, -1.0, 10.0, 20.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, ms, v7, 0.0, 1.0, -1.0, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, ms, v8, 0.0, 1.0, 10.0, -10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, ms, v9, 0.0, 1.0, -100.0, 100.0, 1.0e-3, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, s, v0, 0.0, 1.0, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, s, v1, 0.0, 10.0, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, s, v2, 0.0, 1.0, 10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, s, v3, 0.0, 1.0, 100.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, s, v4, -1.0, 1.0, 10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, s, v5, 0.123456789, 0.987654321, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, s, v6, 1.0, -1.0, 10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, s, v7, 0.0, 1.0, -1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, s, v8, 0.0, 1.0, -10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(trapezoidal, s, v9, 0.0, 1.0, -100.0, 1.0, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, s, v0, 0.0, 1.0, 1.0, 2.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, s, v1, 0.0, 10.0, 1.0, 2.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, s, v2, 0.0, 1.0, 10.0, 20.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, s, v3, 0.0, 1.0, 100.0, 200.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, s, v4, -1.0, 1.0, 10.0, 0.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, s, v5, 0.123456789, 0.987654321, 1.0, 2.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, s, v6, 1.0, -1.0, 10.0, 20.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, s, v7, 0.0, 1.0, -1.0, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, s, v8, 0.0, 1.0, 10.0, -10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(trapezoidal, s, v9, 0.0, 1.0, -100.0, 100.0, 1.0, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, us, v0, 0.0, 1.0, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, us, v1, 0.0, 10.0, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, us, v2, 0.0, 1.0, 10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, us, v3, 0.0, 1.0, 100.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, us, v4, -1.0, 1.0, 10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, us, v5, 0.123456789, 0.987654321, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, us, v6, 1.0, -1.0, 10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, us, v7, 0.0, 1.0, -1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, us, v8, 0.0, 1.0, -10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, us, v9, 0.0, 1.0, -100.0, 1.0e-6, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, us, v0, 0.0, 1.0, 1.0, 2.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, us, v1, 0.0, 10.0, 1.0, 2.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, us, v2, 0.0, 1.0, 10.0, 20.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, us, v3, 0.0, 1.0, 100.0, 200.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, us, v4, -1.0, 1.0, 10.0, 0.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, us, v5, 0.123456789, 0.987654321, 1.0, 2.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, us, v6, 1.0, -1.0, 10.0, 20.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, us, v7, 0.0, 1.0, -1.0, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, us, v8, 0.0, 1.0, 10.0, -10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, us, v9, 0.0, 1.0, -100.0, 100.0, 1.0e-6, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, ms, v0, 0.0, 1.0, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, ms, v1, 0.0, 10.0, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, ms, v2, 0.0, 1.0, 10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, ms, v3, 0.0, 1.0, 100.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, ms, v4, -1.0, 1.0, 10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, ms, v5, 0.123456789, 0.987654321, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, ms, v6, 1.0, -1.0, 10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, ms, v7, 0.0, 1.0, -1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, ms, v8, 0.0, 1.0, -10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, ms, v9, 0.0, 1.0, -100.0, 1.0e-3, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, ms, v0, 0.0, 1.0, 1.0, 2.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, ms, v1, 0.0, 10.0, 1.0, 2.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, ms, v2, 0.0, 1.0, 10.0, 20.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, ms, v3, 0.0, 1.0, 100.0, 200.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, ms, v4, -1.0, 1.0, 10.0, 0.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, ms, v5, 0.123456789, 0.987654321, 1.0, 2.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, ms, v6, 1.0, -1.0, 10.0, 20.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, ms, v7, 0.0, 1.0, -1.0, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, ms, v8, 0.0, 1.0, 10.0, -10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, ms, v9, 0.0, 1.0, -100.0, 100.0, 1.0e-3, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, s, v0, 0.0, 1.0, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, s, v1, 0.0, 10.0, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, s, v2, 0.0, 1.0, 10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, s, v3, 0.0, 1.0, 100.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, s, v4, -1.0, 1.0, 10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, s, v5, 0.123456789, 0.987654321, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, s, v6, 1.0, -1.0, 10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, s, v7, 0.0, 1.0, -1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, s, v8, 0.0, 1.0, -10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons, s, v9, 0.0, 1.0, -100.0, 1.0, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, s, v0, 0.0, 1.0, 1.0, 2.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, s, v1, 0.0, 10.0, 1.0, 2.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, s, v2, 0.0, 1.0, 10.0, 20.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, s, v3, 0.0, 1.0, 100.0, 200.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, s, v4, -1.0, 1.0, 10.0, 0.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, s, v5, 0.123456789, 0.987654321, 1.0, 2.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, s, v6, 1.0, -1.0, 10.0, 20.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, s, v7, 0.0, 1.0, -1.0, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, s, v8, 0.0, 1.0, 10.0, -10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons, s, v9, 0.0, 1.0, -100.0, 100.0, 1.0, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, us, v0, 0.0, 1.0, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, us, v1, 0.0, 10.0, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, us, v2, 0.0, 1.0, 10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, us, v3, 0.0, 1.0, 100.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, us, v4, -1.0, 1.0, 10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, us, v5, 0.123456789, 0.987654321, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, us, v6, 1.0, -1.0, 10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, us, v7, 0.0, 1.0, -1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, us, v8, 0.0, 1.0, -10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, us, v9, 0.0, 1.0, -100.0, 1.0e-6, 1.0e-9, 1.0)

// Note: time versions less accurate
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, us, v0, 0.0, 1.0, 1.0, 2.0, 1.0e-6, 1.0e-9, 1.0e2)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, us, v1, 0.0, 10.0, 1.0, 2.0, 1.0e-6, 1.0e-9, 1.0e2)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, us, v2, 0.0, 1.0, 10.0, 20.0, 1.0e-6, 1.0e-9, 1.0e2)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, us, v3, 0.0, 1.0, 100.0, 200.0, 1.0e-6, 1.0e-9, 1.0e2)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, us, v4, -1.0, 1.0, 10.0, 0.0, 1.0e-6, 1.0e-9, 1.0e2)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, us, v5, 0.123456789, 0.987654321, 1.0, 2.0, 1.0e-6, 1.0e-9, 1.0e2)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, us, v6, 1.0, -1.0, 10.0, 20.0, 1.0e-6, 1.0e-9, 1.0e2)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, us, v7, 0.0, 1.0, -1.0, 1.0, 1.0e-6, 1.0e-9, 1.0e2)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, us, v8, 0.0, 1.0, 10.0, -10.0, 1.0e-6, 1.0e-9, 1.0e2)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, us, v9, 0.0, 1.0, -100.0, 100.0, 1.0e-6, 1.0e-9, 1.0e2)

GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, ms, v0, 0.0, 1.0, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, ms, v1, 0.0, 10.0, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, ms, v2, 0.0, 1.0, 10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, ms, v3, 0.0, 1.0, 100.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, ms, v4, -1.0, 1.0, 10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, ms, v5, 0.123456789, 0.987654321, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, ms, v6, 1.0, -1.0, 10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, ms, v7, 0.0, 1.0, -1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, ms, v8, 0.0, 1.0, -10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, ms, v9, 0.0, 1.0, -100.0, 1.0e-3, 1.0e-9, 1.0)

// Note: time versions less accurate
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, ms, v0, 0.0, 1.0, 1.0, 2.0, 1.0e-3, 1.0e-9, 1.0e2)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, ms, v1, 0.0, 10.0, 1.0, 2.0, 1.0e-3, 1.0e-9, 1.0e2)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, ms, v2, 0.0, 1.0, 10.0, 20.0, 1.0e-3, 1.0e-9, 1.0e2)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, ms, v3, 0.0, 1.0, 100.0, 200.0, 1.0e-3, 1.0e-9, 1.0e2)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, ms, v4, -1.0, 1.0, 10.0, 0.0, 1.0e-3, 1.0e-9, 1.0e2)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, ms, v5, 0.123456789, 0.987654321, 1.0, 2.0, 1.0e-3, 1.0e-9, 1.0e2)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, ms, v6, 1.0, -1.0, 10.0, 20.0, 1.0e-3, 1.0e-9, 1.0e2)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, ms, v7, 0.0, 1.0, -1.0, 1.0, 1.0e-3, 1.0e-9, 1.0e2)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, ms, v8, 0.0, 1.0, 10.0, -10.0, 1.0e-3, 1.0e-9, 1.0e2)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, ms, v9, 0.0, 1.0, -100.0, 100.0, 1.0e-3, 1.0e-9, 1.0e2)

GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, s, v0, 0.0, 1.0, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, s, v1, 0.0, 10.0, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, s, v2, 0.0, 1.0, 10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, s, v3, 0.0, 1.0, 100.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, s, v4, -1.0, 1.0, 10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, s, v5, 0.123456789, 0.987654321, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, s, v6, 1.0, -1.0, 10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, s, v7, 0.0, 1.0, -1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, s, v8, 0.0, 1.0, -10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(simpsons38, s, v9, 0.0, 1.0, -100.0, 1.0, 1.0e-9, 1.0)

// Note: time versions less accurate
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, s, v0, 0.0, 1.0, 1.0, 2.0, 1.0, 1.0e-9, 1.0e2)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, s, v1, 0.0, 10.0, 1.0, 2.0, 1.0, 1.0e-9, 1.0e2)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, s, v2, 0.0, 1.0, 10.0, 20.0, 1.0, 1.0e-9, 1.0e2)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, s, v3, 0.0, 1.0, 100.0, 200.0, 1.0, 1.0e-9, 1.0e2)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, s, v4, -1.0, 1.0, 10.0, 0.0, 1.0, 1.0e-9, 1.0e2)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, s, v5, 0.123456789, 0.987654321, 1.0, 2.0, 1.0, 1.0e-9, 1.0e2)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, s, v6, 1.0, -1.0, 10.0, 20.0, 1.0, 1.0e-9, 1.0e2)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, s, v7, 0.0, 1.0, -1.0, 1.0, 1.0, 1.0e-9, 1.0e2)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, s, v8, 0.0, 1.0, 10.0, -10.0, 1.0, 1.0e-9, 1.0e2)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(simpsons38, s, v9, 0.0, 1.0, -100.0, 100.0, 1.0, 1.0e-9, 1.0e2)

GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, us, v0, 0.0, 1.0, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, us, v1, 0.0, 10.0, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, us, v2, 0.0, 1.0, 10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, us, v3, 0.0, 1.0, 100.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, us, v4, -1.0, 1.0, 10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, us, v5, 0.123456789, 0.987654321, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, us, v6, 1.0, -1.0, 10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, us, v7, 0.0, 1.0, -1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, us, v8, 0.0, 1.0, -10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, us, v9, 0.0, 1.0, -100.0, 1.0e-6, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, us, v0, 0.0, 1.0, 1.0, 2.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, us, v1, 0.0, 10.0, 1.0, 2.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, us, v2, 0.0, 1.0, 10.0, 20.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, us, v3, 0.0, 1.0, 100.0, 200.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, us, v4, -1.0, 1.0, 10.0, 0.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, us, v5, 0.123456789, 0.987654321, 1.0, 2.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, us, v6, 1.0, -1.0, 10.0, 20.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, us, v7, 0.0, 1.0, -1.0, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, us, v8, 0.0, 1.0, 10.0, -10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, us, v9, 0.0, 1.0, -100.0, 100.0, 1.0e-6, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, ms, v0, 0.0, 1.0, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, ms, v1, 0.0, 10.0, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, ms, v2, 0.0, 1.0, 10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, ms, v3, 0.0, 1.0, 100.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, ms, v4, -1.0, 1.0, 10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, ms, v5, 0.123456789, 0.987654321, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, ms, v6, 1.0, -1.0, 10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, ms, v7, 0.0, 1.0, -1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, ms, v8, 0.0, 1.0, -10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, ms, v9, 0.0, 1.0, -100.0, 1.0e-3, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, ms, v0, 0.0, 1.0, 1.0, 2.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, ms, v1, 0.0, 10.0, 1.0, 2.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, ms, v2, 0.0, 1.0, 10.0, 20.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, ms, v3, 0.0, 1.0, 100.0, 200.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, ms, v4, -1.0, 1.0, 10.0, 0.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, ms, v5, 0.123456789, 0.987654321, 1.0, 2.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, ms, v6, 1.0, -1.0, 10.0, 20.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, ms, v7, 0.0, 1.0, -1.0, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, ms, v8, 0.0, 1.0, 10.0, -10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, ms, v9, 0.0, 1.0, -100.0, 100.0, 1.0e-3, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, s, v0, 0.0, 1.0, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, s, v1, 0.0, 10.0, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, s, v2, 0.0, 1.0, 10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, s, v3, 0.0, 1.0, 100.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, s, v4, -1.0, 1.0, 10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, s, v5, 0.123456789, 0.987654321, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, s, v6, 1.0, -1.0, 10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, s, v7, 0.0, 1.0, -1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, s, v8, 0.0, 1.0, -10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(booles, s, v9, 0.0, 1.0, -100.0, 1.0, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, s, v0, 0.0, 1.0, 1.0, 2.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, s, v1, 0.0, 10.0, 1.0, 2.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, s, v2, 0.0, 1.0, 10.0, 20.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, s, v3, 0.0, 1.0, 100.0, 200.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, s, v4, -1.0, 1.0, 10.0, 0.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, s, v5, 0.123456789, 0.987654321, 1.0, 2.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, s, v6, 1.0, -1.0, 10.0, 20.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, s, v7, 0.0, 1.0, -1.0, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, s, v8, 0.0, 1.0, 10.0, -10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(booles, s, v9, 0.0, 1.0, -100.0, 100.0, 1.0, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, us, v0, 0.0, 1.0, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, us, v1, 0.0, 10.0, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, us, v2, 0.0, 1.0, 10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, us, v3, 0.0, 1.0, 100.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, us, v4, -1.0, 1.0, 10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, us, v5, 0.123456789, 0.987654321, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, us, v6, 1.0, -1.0, 10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, us, v7, 0.0, 1.0, -1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, us, v8, 0.0, 1.0, -10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, us, v9, 0.0, 1.0, -100.0, 1.0e-6, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, us, v0, 0.0, 1.0, 1.0, 2.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, us, v1, 0.0, 10.0, 1.0, 2.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, us, v2, 0.0, 1.0, 10.0, 20.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, us, v3, 0.0, 1.0, 100.0, 200.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, us, v4, -1.0, 1.0, 10.0, 0.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, us, v5, 0.123456789, 0.987654321, 1.0, 2.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, us, v6, 1.0, -1.0, 10.0, 20.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, us, v7, 0.0, 1.0, -1.0, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, us, v8, 0.0, 1.0, 10.0, -10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, us, v9, 0.0, 1.0, -100.0, 100.0, 1.0e-6, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, ms, v0, 0.0, 1.0, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, ms, v1, 0.0, 10.0, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, ms, v2, 0.0, 1.0, 10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, ms, v3, 0.0, 1.0, 100.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, ms, v4, -1.0, 1.0, 10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, ms, v5, 0.123456789, 0.987654321, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, ms, v6, 1.0, -1.0, 10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, ms, v7, 0.0, 1.0, -1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, ms, v8, 0.0, 1.0, -10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, ms, v9, 0.0, 1.0, -100.0, 1.0e-3, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, ms, v0, 0.0, 1.0, 1.0, 2.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, ms, v1, 0.0, 10.0, 1.0, 2.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, ms, v2, 0.0, 1.0, 10.0, 20.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, ms, v3, 0.0, 1.0, 100.0, 200.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, ms, v4, -1.0, 1.0, 10.0, 0.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, ms, v5, 0.123456789, 0.987654321, 1.0, 2.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, ms, v6, 1.0, -1.0, 10.0, 20.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, ms, v7, 0.0, 1.0, -1.0, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, ms, v8, 0.0, 1.0, 10.0, -10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, ms, v9, 0.0, 1.0, -100.0, 100.0, 1.0e-3, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, s, v0, 0.0, 1.0, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, s, v1, 0.0, 10.0, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, s, v2, 0.0, 1.0, 10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, s, v3, 0.0, 1.0, 100.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, s, v4, -1.0, 1.0, 10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, s, v5, 0.123456789, 0.987654321, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, s, v6, 1.0, -1.0, 10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, s, v7, 0.0, 1.0, -1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, s, v8, 0.0, 1.0, -10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(rectangle, s, v9, 0.0, 1.0, -100.0, 1.0, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, s, v0, 0.0, 1.0, 1.0, 2.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, s, v1, 0.0, 10.0, 1.0, 2.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, s, v2, 0.0, 1.0, 10.0, 20.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, s, v3, 0.0, 1.0, 100.0, 200.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, s, v4, -1.0, 1.0, 10.0, 0.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, s, v5, 0.123456789, 0.987654321, 1.0, 2.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, s, v6, 1.0, -1.0, 10.0, 20.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, s, v7, 0.0, 1.0, -1.0, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, s, v8, 0.0, 1.0, 10.0, -10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(rectangle, s, v9, 0.0, 1.0, -100.0, 100.0, 1.0, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, us, v0, 0.0, 1.0, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, us, v1, 0.0, 10.0, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, us, v2, 0.0, 1.0, 10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, us, v3, 0.0, 1.0, 100.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, us, v4, -1.0, 1.0, 10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, us, v5, 0.123456789, 0.987654321, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, us, v6, 1.0, -1.0, 10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, us, v7, 0.0, 1.0, -1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, us, v8, 0.0, 1.0, -10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, us, v9, 0.0, 1.0, -100.0, 1.0e-6, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, us, v0, 0.0, 1.0, 1.0, 2.0, 1.0e-6, 1.0e-9, 1.0e1)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, us, v1, 0.0, 10.0, 1.0, 2.0, 1.0e-6, 1.0e-9, 1.0e1)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, us, v2, 0.0, 1.0, 10.0, 20.0, 1.0e-6, 1.0e-9, 1.0e1)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, us, v3, 0.0, 1.0, 100.0, 200.0, 1.0e-6, 1.0e-9, 1.0e1)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, us, v4, -1.0, 1.0, 10.0, 0.0, 1.0e-6, 1.0e-9, 1.0e1)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, us, v5, 0.123456789, 0.987654321, 1.0, 2.0, 1.0e-6, 1.0e-9, 1.0e1)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, us, v6, 1.0, -1.0, 10.0, 20.0, 1.0e-6, 1.0e-9, 1.0e1)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, us, v7, 0.0, 1.0, -1.0, 1.0, 1.0e-6, 1.0e-9, 1.0e1)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, us, v8, 0.0, 1.0, 10.0, -10.0, 1.0e-6, 1.0e-9, 1.0e1)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, us, v9, 0.0, 1.0, -100.0, 100.0, 1.0e-6, 1.0e-9, 1.0e1)

GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, ms, v0, 0.0, 1.0, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, ms, v1, 0.0, 10.0, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, ms, v2, 0.0, 1.0, 10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, ms, v3, 0.0, 1.0, 100.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, ms, v4, -1.0, 1.0, 10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, ms, v5, 0.123456789, 0.987654321, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, ms, v6, 1.0, -1.0, 10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, ms, v7, 0.0, 1.0, -1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, ms, v8, 0.0, 1.0, -10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, ms, v9, 0.0, 1.0, -100.0, 1.0e-3, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, ms, v0, 0.0, 1.0, 1.0, 2.0, 1.0e-3, 1.0e-9, 1.0e1)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, ms, v1, 0.0, 10.0, 1.0, 2.0, 1.0e-3, 1.0e-9, 1.0e1)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, ms, v2, 0.0, 1.0, 10.0, 20.0, 1.0e-3, 1.0e-9, 1.0e1)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, ms, v3, 0.0, 1.0, 100.0, 200.0, 1.0e-3, 1.0e-9, 1.0e1)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, ms, v4, -1.0, 1.0, 10.0, 0.0, 1.0e-3, 1.0e-9, 1.0e1)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, ms, v5, 0.123456789, 0.987654321, 1.0, 2.0, 1.0e-3, 1.0e-9, 1.0e1)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, ms, v6, 1.0, -1.0, 10.0, 20.0, 1.0e-3, 1.0e-9, 1.0e1)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, ms, v7, 0.0, 1.0, -1.0, 1.0, 1.0e-3, 1.0e-9, 1.0e1)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, ms, v8, 0.0, 1.0, 10.0, -10.0, 1.0e-3, 1.0e-9, 1.0e1)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, ms, v9, 0.0, 1.0, -100.0, 100.0, 1.0e-3, 1.0e-9, 1.0e1)

GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, s, v0, 0.0, 1.0, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, s, v1, 0.0, 10.0, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, s, v2, 0.0, 1.0, 10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, s, v3, 0.0, 1.0, 100.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, s, v4, -1.0, 1.0, 10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, s, v5, 0.123456789, 0.987654321, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, s, v6, 1.0, -1.0, 10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, s, v7, 0.0, 1.0, -1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, s, v8, 0.0, 1.0, -10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open1, s, v9, 0.0, 1.0, -100.0, 1.0, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, s, v0, 0.0, 1.0, 1.0, 2.0, 1.0, 1.0e-9, 1.0e1)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, s, v1, 0.0, 10.0, 1.0, 2.0, 1.0, 1.0e-9, 1.0e1)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, s, v2, 0.0, 1.0, 10.0, 20.0, 1.0, 1.0e-9, 1.0e1)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, s, v3, 0.0, 1.0, 100.0, 200.0, 1.0, 1.0e-9, 1.0e1)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, s, v4, -1.0, 1.0, 10.0, 0.0, 1.0, 1.0e-9, 1.0e1)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, s, v5, 0.123456789, 0.987654321, 1.0, 2.0, 1.0, 1.0e-9, 1.0e1)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, s, v6, 1.0, -1.0, 10.0, 20.0, 1.0, 1.0e-9, 1.0e1)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, s, v7, 0.0, 1.0, -1.0, 1.0, 1.0, 1.0e-9, 1.0e1)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, s, v8, 0.0, 1.0, 10.0, -10.0, 1.0, 1.0e-9, 1.0e1)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open1, s, v9, 0.0, 1.0, -100.0, 100.0, 1.0, 1.0e-9, 1.0e1)

GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, us, v0, 0.0, 1.0, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, us, v1, 0.0, 10.0, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, us, v2, 0.0, 1.0, 10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, us, v3, 0.0, 1.0, 100.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, us, v4, -1.0, 1.0, 10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, us, v5, 0.123456789, 0.987654321, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, us, v6, 1.0, -1.0, 10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, us, v7, 0.0, 1.0, -1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, us, v8, 0.0, 1.0, -10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, us, v9, 0.0, 1.0, -100.0, 1.0e-6, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, us, v0, 0.0, 1.0, 1.0, 2.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, us, v1, 0.0, 10.0, 1.0, 2.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, us, v2, 0.0, 1.0, 10.0, 20.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, us, v3, 0.0, 1.0, 100.0, 200.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, us, v4, -1.0, 1.0, 10.0, 0.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, us, v5, 0.123456789, 0.987654321, 1.0, 2.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, us, v6, 1.0, -1.0, 10.0, 20.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, us, v7, 0.0, 1.0, -1.0, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, us, v8, 0.0, 1.0, 10.0, -10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, us, v9, 0.0, 1.0, -100.0, 100.0, 1.0e-6, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, ms, v0, 0.0, 1.0, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, ms, v1, 0.0, 10.0, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, ms, v2, 0.0, 1.0, 10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, ms, v3, 0.0, 1.0, 100.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, ms, v4, -1.0, 1.0, 10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, ms, v5, 0.123456789, 0.987654321, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, ms, v6, 1.0, -1.0, 10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, ms, v7, 0.0, 1.0, -1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, ms, v8, 0.0, 1.0, -10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, ms, v9, 0.0, 1.0, -100.0, 1.0e-3, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, ms, v0, 0.0, 1.0, 1.0, 2.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, ms, v1, 0.0, 10.0, 1.0, 2.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, ms, v2, 0.0, 1.0, 10.0, 20.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, ms, v3, 0.0, 1.0, 100.0, 200.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, ms, v4, -1.0, 1.0, 10.0, 0.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, ms, v5, 0.123456789, 0.987654321, 1.0, 2.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, ms, v6, 1.0, -1.0, 10.0, 20.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, ms, v7, 0.0, 1.0, -1.0, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, ms, v8, 0.0, 1.0, 10.0, -10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, ms, v9, 0.0, 1.0, -100.0, 100.0, 1.0e-3, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, s, v0, 0.0, 1.0, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, s, v1, 0.0, 10.0, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, s, v2, 0.0, 1.0, 10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, s, v3, 0.0, 1.0, 100.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, s, v4, -1.0, 1.0, 10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, s, v5, 0.123456789, 0.987654321, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, s, v6, 1.0, -1.0, 10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, s, v7, 0.0, 1.0, -1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, s, v8, 0.0, 1.0, -10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(milnes, s, v9, 0.0, 1.0, -100.0, 1.0, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, s, v0, 0.0, 1.0, 1.0, 2.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, s, v1, 0.0, 10.0, 1.0, 2.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, s, v2, 0.0, 1.0, 10.0, 20.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, s, v3, 0.0, 1.0, 100.0, 200.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, s, v4, -1.0, 1.0, 10.0, 0.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, s, v5, 0.123456789, 0.987654321, 1.0, 2.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, s, v6, 1.0, -1.0, 10.0, 20.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, s, v7, 0.0, 1.0, -1.0, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, s, v8, 0.0, 1.0, 10.0, -10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(milnes, s, v9, 0.0, 1.0, -100.0, 100.0, 1.0, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, us, v0, 0.0, 1.0, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, us, v1, 0.0, 10.0, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, us, v2, 0.0, 1.0, 10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, us, v3, 0.0, 1.0, 100.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, us, v4, -1.0, 1.0, 10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, us, v5, 0.123456789, 0.987654321, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, us, v6, 1.0, -1.0, 10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, us, v7, 0.0, 1.0, -1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, us, v8, 0.0, 1.0, -10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, us, v9, 0.0, 1.0, -100.0, 1.0e-6, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, us, v0, 0.0, 1.0, 1.0, 2.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, us, v1, 0.0, 10.0, 1.0, 2.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, us, v2, 0.0, 1.0, 10.0, 20.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, us, v3, 0.0, 1.0, 100.0, 200.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, us, v4, -1.0, 1.0, 10.0, 0.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, us, v5, 0.123456789, 0.987654321, 1.0, 2.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, us, v6, 1.0, -1.0, 10.0, 20.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, us, v7, 0.0, 1.0, -1.0, 1.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, us, v8, 0.0, 1.0, 10.0, -10.0, 1.0e-6, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, us, v9, 0.0, 1.0, -100.0, 100.0, 1.0e-6, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, ms, v0, 0.0, 1.0, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, ms, v1, 0.0, 10.0, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, ms, v2, 0.0, 1.0, 10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, ms, v3, 0.0, 1.0, 100.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, ms, v4, -1.0, 1.0, 10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, ms, v5, 0.123456789, 0.987654321, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, ms, v6, 1.0, -1.0, 10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, ms, v7, 0.0, 1.0, -1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, ms, v8, 0.0, 1.0, -10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, ms, v9, 0.0, 1.0, -100.0, 1.0e-3, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, ms, v0, 0.0, 1.0, 1.0, 2.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, ms, v1, 0.0, 10.0, 1.0, 2.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, ms, v2, 0.0, 1.0, 10.0, 20.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, ms, v3, 0.0, 1.0, 100.0, 200.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, ms, v4, -1.0, 1.0, 10.0, 0.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, ms, v5, 0.123456789, 0.987654321, 1.0, 2.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, ms, v6, 1.0, -1.0, 10.0, 20.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, ms, v7, 0.0, 1.0, -1.0, 1.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, ms, v8, 0.0, 1.0, 10.0, -10.0, 1.0e-3, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, ms, v9, 0.0, 1.0, -100.0, 100.0, 1.0e-3, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, s, v0, 0.0, 1.0, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, s, v1, 0.0, 10.0, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, s, v2, 0.0, 1.0, 10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, s, v3, 0.0, 1.0, 100.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, s, v4, -1.0, 1.0, 10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, s, v5, 0.123456789, 0.987654321, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, s, v6, 1.0, -1.0, 10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, s, v7, 0.0, 1.0, -1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, s, v8, 0.0, 1.0, -10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_CONSTANT_STEP(open3, s, v9, 0.0, 1.0, -100.0, 1.0, 1.0e-9, 1.0)

GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, s, v0, 0.0, 1.0, 1.0, 2.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, s, v1, 0.0, 10.0, 1.0, 2.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, s, v2, 0.0, 1.0, 10.0, 20.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, s, v3, 0.0, 1.0, 100.0, 200.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, s, v4, -1.0, 1.0, 10.0, 0.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, s, v5, 0.123456789, 0.987654321, 1.0, 2.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, s, v6, 1.0, -1.0, 10.0, 20.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, s, v7, 0.0, 1.0, -1.0, 1.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, s, v8, 0.0, 1.0, 10.0, -10.0, 1.0, 1.0e-9, 1.0)
GENERATE_INTEGRATION_TEST_LINEAR_STEP(open3, s, v9, 0.0, 1.0, -100.0, 100.0, 1.0, 1.0e-9, 1.0)
