#include "mathbox/integrator.hpp"

#include <gtest/gtest.h>

template<typename Integrator>
void test_integrate_constant(const int num_subintervals) {
    Integrator integrator;
    const std::vector<double> starts = {{0.0, -1.0, 0.9}};
    const std::vector<double> ends = {{1.0, 1.4, 7.7}};
    const std::vector<double> constants = {{3.0, -1.7, 0.0}};
    const std::vector<math::IntegrandFunction<double>> integrands = {{+[](const double) -> double { return 0.0; },
            +[](const double) -> double { return 1.0; }, +[](const double) -> double { return 4.3; }}};
    for (const double start : starts) {
        for (const double end : ends) {
            for (const auto integrand : integrands) {
                EXPECT_NEAR(integrator.integrate(start, end, num_subintervals, integrand),
                        (end - start) * integrand(start), 1.0e-12);
            }
            for (const auto constant : constants) {
                EXPECT_NEAR(integrator.integrate(start, end, num_subintervals,
                                    [constant](const double) -> double { return constant; }),
                        (end - start) * constant, 1.0e-12);
            }
        }
    }
}

template<typename Integrator>
void test_integrate_linear(const int num_subintervals) {
    Integrator integrator;
    const std::vector<double> starts = {{0.0, -1.0, 0.9}};
    const std::vector<double> ends = {{1.0, 1.4, 7.7}};
    const std::vector<double> start_heights = {{0.0, 1.0, 4.3}};
    const std::vector<double> end_heights = {{-0.5, 3.1, 9.7}};
    for (const double start : starts) {
        for (const double end : ends) {
            for (const double start_height : start_heights) {
                for (const double end_height : end_heights) {
                    EXPECT_NEAR(integrator.integrate(start, end, num_subintervals,
                                        std::bind(math::linear_function<double>, std::placeholders::_1, start, end,
                                                start_height, end_height)),
                            (end - start) * (start_height + end_height) / 2.0, 1.0e-12);
                    EXPECT_NEAR(integrator.integrate(start, end, num_subintervals,
                                        [start, end, start_height, end_height](const double t) -> double {
                                            return math::linear_function(t, start, end, start_height, end_height);
                                        }),
                            (end - start) * (start_height + end_height) / 2.0, 1.0e-12);
                }
            }
        }
    }
}

TEST(TrapezoidalIntegration, integrate_constant) {
    test_integrate_constant<math::newton_cotes::trapezoidal::Integrator<double>>(1);
    test_integrate_constant<math::newton_cotes::trapezoidal::Integrator<double>>(13);
    test_integrate_constant<math::newton_cotes::trapezoidal::Integrator<double>>(100);
}

TEST(TrapezoidalIntegration, integrate_linear) {
    test_integrate_linear<math::newton_cotes::trapezoidal::Integrator<double>>(1);
    test_integrate_linear<math::newton_cotes::trapezoidal::Integrator<double>>(13);
    test_integrate_linear<math::newton_cotes::trapezoidal::Integrator<double>>(100);
}

TEST(SimpsonIntegration, integrate_constant) {
    test_integrate_constant<math::newton_cotes::simpsons::Integrator<double>>(1);
    test_integrate_constant<math::newton_cotes::simpsons::Integrator<double>>(13);
    test_integrate_constant<math::newton_cotes::simpsons::Integrator<double>>(100);
}

TEST(SimpsonIntegration, integrate_linear) {
    test_integrate_linear<math::newton_cotes::simpsons::Integrator<double>>(1);
    test_integrate_linear<math::newton_cotes::simpsons::Integrator<double>>(13);
    test_integrate_linear<math::newton_cotes::simpsons::Integrator<double>>(100);
}

TEST(Simpsons38Integration, integrate_constant) {
    test_integrate_constant<math::newton_cotes::simpsons38::Integrator<double>>(1);
    test_integrate_constant<math::newton_cotes::simpsons38::Integrator<double>>(13);
    test_integrate_constant<math::newton_cotes::simpsons38::Integrator<double>>(100);
}

TEST(Simpsons38Integration, integrate_linear) {
    test_integrate_linear<math::newton_cotes::simpsons38::Integrator<double>>(1);
    test_integrate_linear<math::newton_cotes::simpsons38::Integrator<double>>(13);
    test_integrate_linear<math::newton_cotes::simpsons38::Integrator<double>>(100);
}

TEST(BoolesIntegration, integrate_constant) {
    test_integrate_constant<math::newton_cotes::booles::Integrator<double>>(1);
    test_integrate_constant<math::newton_cotes::booles::Integrator<double>>(13);
    test_integrate_constant<math::newton_cotes::booles::Integrator<double>>(100);
}

TEST(BoolesIntegration, integrate_linear) {
    test_integrate_linear<math::newton_cotes::booles::Integrator<double>>(1);
    test_integrate_linear<math::newton_cotes::booles::Integrator<double>>(13);
    test_integrate_linear<math::newton_cotes::booles::Integrator<double>>(100);
}

TEST(RectangleIntegration, integrate_constant) {
    test_integrate_constant<math::newton_cotes::rectangle::Integrator<double>>(1);
    test_integrate_constant<math::newton_cotes::rectangle::Integrator<double>>(13);
    test_integrate_constant<math::newton_cotes::rectangle::Integrator<double>>(100);
}

TEST(RectangleIntegration, integrate_linear) {
    test_integrate_linear<math::newton_cotes::rectangle::Integrator<double>>(1);
    test_integrate_linear<math::newton_cotes::rectangle::Integrator<double>>(13);
    test_integrate_linear<math::newton_cotes::rectangle::Integrator<double>>(100);
}

TEST(Open1Integration, integrate_constant) {
    test_integrate_constant<math::newton_cotes::open1::Integrator<double>>(1);
    test_integrate_constant<math::newton_cotes::open1::Integrator<double>>(13);
    test_integrate_constant<math::newton_cotes::open1::Integrator<double>>(100);
}

TEST(Open1Integration, integrate_linear) {
    test_integrate_linear<math::newton_cotes::open1::Integrator<double>>(1);
    test_integrate_linear<math::newton_cotes::open1::Integrator<double>>(13);
    test_integrate_linear<math::newton_cotes::open1::Integrator<double>>(100);
}

TEST(MilnesIntegration, integrate_constant) {
    test_integrate_constant<math::newton_cotes::milnes::Integrator<double>>(1);
    test_integrate_constant<math::newton_cotes::milnes::Integrator<double>>(13);
    test_integrate_constant<math::newton_cotes::milnes::Integrator<double>>(100);
}

TEST(MilnesIntegration, integrate_linear) {
    test_integrate_linear<math::newton_cotes::milnes::Integrator<double>>(1);
    test_integrate_linear<math::newton_cotes::milnes::Integrator<double>>(13);
    test_integrate_linear<math::newton_cotes::milnes::Integrator<double>>(100);
}

TEST(Open3Integration, integrate_constant) {
    test_integrate_constant<math::newton_cotes::open3::Integrator<double>>(1);
    test_integrate_constant<math::newton_cotes::open3::Integrator<double>>(13);
    test_integrate_constant<math::newton_cotes::open3::Integrator<double>>(100);
}

TEST(Open3Integration, integrate_linear) {
    test_integrate_linear<math::newton_cotes::open3::Integrator<double>>(1);
    test_integrate_linear<math::newton_cotes::open3::Integrator<double>>(13);
    test_integrate_linear<math::newton_cotes::open3::Integrator<double>>(100);
}
