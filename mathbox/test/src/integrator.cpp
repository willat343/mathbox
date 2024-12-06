#include "mathbox/integrator.hpp"

#include <gtest/gtest.h>

#include <chrono>

template<typename Integrator>
void test_integrate_constant(const std::vector<typename Integrator::IndependentVariableType>& starts,
        const std::vector<typename Integrator::IndependentVariableType>& ends,
        const std::vector<typename Integrator::ArithmeticType> constants, const int num_subintervals,
        const double precision = 1.0e-12) {
    Integrator integrator;
    using ArithmeticType = Integrator::ArithmeticType;
    using IndependentVariableType = Integrator::IndependentVariableType;
    for (const IndependentVariableType start : starts) {
        for (const IndependentVariableType end : ends) {
            for (const ArithmeticType constant : constants) {
                ArithmeticType expected_value = constant;
                if constexpr (math::is_time_point_v<IndependentVariableType>) {
                    expected_value *= math::to_sec(end - start);
                } else {
                    expected_value *= (end - start);
                }
                EXPECT_NEAR(integrator.integrate(start, end, num_subintervals,
                                    [constant](const IndependentVariableType) -> ArithmeticType { return constant; }),
                        expected_value, precision);
            }
        }
    }
}

template<typename Integrator>
void test_integrate_linear(const std::vector<typename Integrator::IndependentVariableType>& starts,
        const std::vector<typename Integrator::IndependentVariableType>& ends,
        const std::vector<typename Integrator::ArithmeticType> start_heights,
        const std::vector<typename Integrator::ArithmeticType> end_heights, const int num_subintervals,
        const double precision = 1.0e-12) {
    Integrator integrator;
    using ArithmeticType = Integrator::ArithmeticType;
    using IndependentVariableType = Integrator::IndependentVariableType;
    for (const IndependentVariableType start : starts) {
        for (const IndependentVariableType end : ends) {
            for (const ArithmeticType start_height : start_heights) {
                for (const ArithmeticType end_height : end_heights) {
                    typename Integrator::ArithmeticType expected_value = (start_height + end_height) / 2.0;
                    if constexpr (math::is_time_point_v<typename Integrator::IndependentVariableType>) {
                        expected_value *= math::to_sec(end - start);
                    } else {
                        expected_value *= (end - start);
                    }
                    EXPECT_NEAR(integrator.integrate(start, end, num_subintervals,
                                        std::bind(math::linear_function<ArithmeticType, IndependentVariableType>,
                                                std::placeholders::_1, start, end, start_height, end_height)),
                            expected_value, precision);
                    EXPECT_NEAR(integrator.integrate(start, end, num_subintervals,
                                        [start, end, start_height, end_height](
                                                const IndependentVariableType t) -> ArithmeticType {
                                            return math::linear_function(t, start, end, start_height, end_height);
                                        }),
                            expected_value, precision);
                }
            }
        }
    }
}

template<typename IndependentVariableType>
std::vector<IndependentVariableType> start_instances() {
    static_assert("IndependentVariableType not one of specialised types.");
}

template<>
std::vector<double> start_instances<double>() {
    return {{0.0, -1.0, 0.9}};
}

template<>
std::vector<std::chrono::time_point<std::chrono::steady_clock>>
start_instances<std::chrono::time_point<std::chrono::steady_clock>>() {
    return {{std::chrono::time_point<std::chrono::steady_clock>(std::chrono::steady_clock::duration(0)),
            std::chrono::time_point<std::chrono::steady_clock>(std::chrono::steady_clock::duration(10)),
            std::chrono::time_point<std::chrono::steady_clock>(std::chrono::steady_clock::duration(123))}};
}

template<typename IndependentVariableType>
std::vector<IndependentVariableType> end_instances() {
    static_assert("IndependentVariableType not one of specialised types.");
}

template<>
std::vector<double> end_instances<double>() {
    return {{1.0, 1.4, 7.7}};
}

template<>
std::vector<std::chrono::time_point<std::chrono::steady_clock>>
end_instances<std::chrono::time_point<std::chrono::steady_clock>>() {
    return {{std::chrono::time_point<std::chrono::steady_clock>(std::chrono::steady_clock::duration(200)),
            std::chrono::time_point<std::chrono::steady_clock>(std::chrono::steady_clock::duration(1234)),
            std::chrono::time_point<std::chrono::steady_clock>(std::chrono::steady_clock::duration(94357128))}};
}

template<typename ArithmeticType>
std::vector<ArithmeticType> constants_instances(const unsigned int i = 0) {
    static_assert("ArithmeticType not one of specialised types.");
}

template<>
std::vector<double> constants_instances<double>(const unsigned int i) {
    if (i == 0) {
        return {{3.0, -1.7, 0.0}};
    } else {
        return {{-0.5, 3.1, 9.7}};
    }
}

template<typename ArithmeticType, typename IndependentVariableType>
std::vector<math::IntegrandFunction<ArithmeticType, IndependentVariableType>> integrands_instances() {
    static_assert("IntegrandFunction not one of specialised types.");
}

template<typename Integrator>
void test_integrate_constant_with_instances(const double precision = 1.0e-12) {
    using ArithmeticType = typename Integrator::ArithmeticType;
    using IndependentVariableType = typename Integrator::IndependentVariableType;
    const std::vector<IndependentVariableType> starts = start_instances<IndependentVariableType>();
    const std::vector<IndependentVariableType> ends = end_instances<IndependentVariableType>();
    const std::vector<ArithmeticType> constants = constants_instances<ArithmeticType>();
    test_integrate_constant<Integrator>(starts, ends, constants, 1, precision);
    test_integrate_constant<Integrator>(starts, ends, constants, 13, precision);
    test_integrate_constant<Integrator>(starts, ends, constants, 77, precision);
    test_integrate_constant<Integrator>(starts, ends, constants, 100, precision);
    test_integrate_constant<Integrator>(starts, ends, constants, 1000, precision);
}

template<typename Integrator>
void test_integrate_linear_with_instances(const double precision = 1.0e-12) {
    using ArithmeticType = typename Integrator::ArithmeticType;
    using IndependentVariableType = typename Integrator::IndependentVariableType;
    const std::vector<IndependentVariableType> starts = start_instances<IndependentVariableType>();
    const std::vector<IndependentVariableType> ends = end_instances<IndependentVariableType>();
    const std::vector<ArithmeticType> start_heights = constants_instances<ArithmeticType>(0);
    const std::vector<ArithmeticType> end_heights = constants_instances<ArithmeticType>(1);
    test_integrate_linear<Integrator>(starts, ends, start_heights, end_heights, 1, precision);
    test_integrate_linear<Integrator>(starts, ends, start_heights, end_heights, 13, precision);
    test_integrate_linear<Integrator>(starts, ends, start_heights, end_heights, 77, precision);
    test_integrate_linear<Integrator>(starts, ends, start_heights, end_heights, 100, precision);
    test_integrate_linear<Integrator>(starts, ends, start_heights, end_heights, 1000, precision);
}

TEST(TrapezoidalIntegration, integrate_constant) {
    test_integrate_constant_with_instances<math::newton_cotes::trapezoidal::Integrator<double, double>>();
    test_integrate_constant_with_instances<
            math::newton_cotes::trapezoidal::Integrator<double, std::chrono::time_point<std::chrono::steady_clock>>>(
            1.0e-8);
}

TEST(TrapezoidalIntegration, integrate_linear) {
    test_integrate_linear_with_instances<math::newton_cotes::trapezoidal::Integrator<double, double>>();
    test_integrate_linear_with_instances<
            math::newton_cotes::trapezoidal::Integrator<double, std::chrono::time_point<std::chrono::steady_clock>>>(
            1.0e-8);
}

TEST(SimpsonIntegration, integrate_constant) {
    test_integrate_constant_with_instances<math::newton_cotes::simpsons::Integrator<double, double>>();
    test_integrate_constant_with_instances<
            math::newton_cotes::simpsons::Integrator<double, std::chrono::time_point<std::chrono::steady_clock>>>(
            1.0e-8);
}

TEST(SimpsonIntegration, integrate_linear) {
    test_integrate_linear_with_instances<math::newton_cotes::simpsons::Integrator<double, double>>();
    test_integrate_linear_with_instances<
            math::newton_cotes::simpsons::Integrator<double, std::chrono::time_point<std::chrono::steady_clock>>>(
            1.0e-8);
}

TEST(Simpsons38Integration, integrate_constant) {
    test_integrate_constant_with_instances<math::newton_cotes::simpsons38::Integrator<double, double>>();
    test_integrate_constant_with_instances<
            math::newton_cotes::simpsons38::Integrator<double, std::chrono::time_point<std::chrono::steady_clock>>>(
            1.0e-8);
}

TEST(Simpsons38Integration, integrate_linear) {
    test_integrate_linear_with_instances<math::newton_cotes::simpsons38::Integrator<double, double>>();
    test_integrate_linear_with_instances<
            math::newton_cotes::simpsons38::Integrator<double, std::chrono::time_point<std::chrono::steady_clock>>>(
            1.0e-8);
}

TEST(BoolesIntegration, integrate_constant) {
    test_integrate_constant_with_instances<math::newton_cotes::booles::Integrator<double, double>>();
    test_integrate_constant_with_instances<
            math::newton_cotes::booles::Integrator<double, std::chrono::time_point<std::chrono::steady_clock>>>(1.0e-8);
}

TEST(BoolesIntegration, integrate_linear) {
    test_integrate_linear_with_instances<math::newton_cotes::booles::Integrator<double, double>>();
    test_integrate_linear_with_instances<
            math::newton_cotes::booles::Integrator<double, std::chrono::time_point<std::chrono::steady_clock>>>(1.0e-8);
}

TEST(RectangleIntegration, integrate_constant) {
    test_integrate_constant_with_instances<math::newton_cotes::rectangle::Integrator<double, double>>();
    test_integrate_constant_with_instances<
            math::newton_cotes::rectangle::Integrator<double, std::chrono::time_point<std::chrono::steady_clock>>>(
            1.0e-8);
}

TEST(RectangleIntegration, integrate_linear) {
    test_integrate_linear_with_instances<math::newton_cotes::rectangle::Integrator<double, double>>();
    test_integrate_linear_with_instances<
            math::newton_cotes::rectangle::Integrator<double, std::chrono::time_point<std::chrono::steady_clock>>>(
            1.0e-8);
}

TEST(Open1Integration, integrate_constant) {
    test_integrate_constant_with_instances<math::newton_cotes::open1::Integrator<double, double>>();
    test_integrate_constant_with_instances<
            math::newton_cotes::open1::Integrator<double, std::chrono::time_point<std::chrono::steady_clock>>>(1.0e-8);
}

TEST(Open1Integration, integrate_linear) {
    test_integrate_linear_with_instances<math::newton_cotes::open1::Integrator<double, double>>();
    test_integrate_linear_with_instances<
            math::newton_cotes::open1::Integrator<double, std::chrono::time_point<std::chrono::steady_clock>>>(1.0e-8);
}

TEST(MilnesIntegration, integrate_constant) {
    test_integrate_constant_with_instances<math::newton_cotes::milnes::Integrator<double, double>>();
    test_integrate_constant_with_instances<
            math::newton_cotes::milnes::Integrator<double, std::chrono::time_point<std::chrono::steady_clock>>>(1.0e-8);
}

TEST(MilnesIntegration, integrate_linear) {
    test_integrate_linear_with_instances<math::newton_cotes::milnes::Integrator<double, double>>();
    test_integrate_linear_with_instances<
            math::newton_cotes::milnes::Integrator<double, std::chrono::time_point<std::chrono::steady_clock>>>(1.0e-8);
}

TEST(Open3Integration, integrate_constant) {
    test_integrate_constant_with_instances<math::newton_cotes::open3::Integrator<double, double>>();
    test_integrate_constant_with_instances<
            math::newton_cotes::open3::Integrator<double, std::chrono::time_point<std::chrono::steady_clock>>>(1.0e-8);
}

TEST(Open3Integration, integrate_linear) {
    test_integrate_linear_with_instances<math::newton_cotes::open3::Integrator<double, double>>();
    test_integrate_linear_with_instances<
            math::newton_cotes::open3::Integrator<double, std::chrono::time_point<std::chrono::steady_clock>>>(1.0e-8);
}
