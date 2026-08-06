#include "mathbox/statistics.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <vector>

TEST(Statistics, default_construction) {
    const math::Statistics<double> statistics;
    EXPECT_DOUBLE_EQ(statistics.mean(), 0.0);
    EXPECT_DOUBLE_EQ(statistics.population_variance(), 0.0);
    EXPECT_DOUBLE_EQ(statistics.minimum(), 0.0);
    EXPECT_DOUBLE_EQ(statistics.maximum(), 0.0);
}

TEST(Statistics, mean_and_variance_construction) {
    const math::Statistics<double> statistics(5.0, 4.0);
    EXPECT_DOUBLE_EQ(statistics.mean(), 5.0);
    EXPECT_DOUBLE_EQ(statistics.population_variance(), 4.0);
    EXPECT_DOUBLE_EQ(statistics.minimum(), 5.0);
    EXPECT_DOUBLE_EQ(statistics.maximum(), 5.0);
    EXPECT_DOUBLE_EQ(statistics.population_stddev(), 2.0);
}

TEST(Statistics, full_construction) {
    const math::Statistics<double> statistics(5.0, 4.0, 2.0, 9.0);
    EXPECT_DOUBLE_EQ(statistics.mean(), 5.0);
    EXPECT_DOUBLE_EQ(statistics.population_variance(), 4.0);
    EXPECT_DOUBLE_EQ(statistics.minimum(), 2.0);
    EXPECT_DOUBLE_EQ(statistics.maximum(), 9.0);
}

TEST(Statistics, mutable_accessors) {
    math::Statistics<double> statistics;
    statistics.mean() = 1.0;
    statistics.population_variance() = 2.0;
    statistics.minimum() = -1.0;
    statistics.maximum() = 3.0;
    EXPECT_DOUBLE_EQ(statistics.mean(), 1.0);
    EXPECT_DOUBLE_EQ(statistics.population_variance(), 2.0);
    EXPECT_DOUBLE_EQ(statistics.minimum(), -1.0);
    EXPECT_DOUBLE_EQ(statistics.maximum(), 3.0);
}

TEST(RunningStatistics, default_construction) {
    const math::RunningStatistics<double> statistics;
    EXPECT_EQ(statistics.num_samples(), 0u);
    EXPECT_DOUBLE_EQ(statistics.sum(), 0.0);
    EXPECT_DOUBLE_EQ(statistics.sum_of_squares(), 0.0);
    EXPECT_DOUBLE_EQ(statistics.sum_of_square_differences(), 0.0);
    EXPECT_DOUBLE_EQ(statistics.rms(), 0.0);
}

TEST(RunningStatistics, single_sample_construction) {
    const math::RunningStatistics<double> statistics(5.0);
    EXPECT_EQ(statistics.num_samples(), 1u);
    EXPECT_DOUBLE_EQ(statistics.mean(), 5.0);
    EXPECT_DOUBLE_EQ(statistics.population_variance(), 0.0);
    EXPECT_DOUBLE_EQ(statistics.minimum(), 5.0);
    EXPECT_DOUBLE_EQ(statistics.maximum(), 5.0);
    EXPECT_DOUBLE_EQ(statistics.sum(), 5.0);
    EXPECT_DOUBLE_EQ(statistics.sum_of_squares(), 25.0);
    EXPECT_DOUBLE_EQ(statistics.rms(), 5.0);
    EXPECT_DOUBLE_EQ(statistics.sample_variance(), 0.0);
    EXPECT_DOUBLE_EQ(statistics.sample_stddev(), 0.0);
}

TEST(RunningStatistics, update_matches_known_dataset) {
    // Classic textbook dataset: {2, 4, 4, 4, 5, 5, 7, 9}, mean = 5, population variance = 4.
    math::RunningStatistics<double> statistics;
    for (const double sample : {2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0}) {
        statistics.update(sample);
    }
    EXPECT_EQ(statistics.num_samples(), 8u);
    EXPECT_DOUBLE_EQ(statistics.mean(), 5.0);
    EXPECT_DOUBLE_EQ(statistics.minimum(), 2.0);
    EXPECT_DOUBLE_EQ(statistics.maximum(), 9.0);
    EXPECT_DOUBLE_EQ(statistics.sum(), 40.0);
    EXPECT_NEAR(statistics.sum_of_square_differences(), 32.0, 1.0e-9);
    EXPECT_NEAR(statistics.sum_of_squares(), 232.0, 1.0e-9);
    EXPECT_NEAR(statistics.population_variance(), 4.0, 1.0e-9);
    EXPECT_NEAR(statistics.sample_variance(), 32.0 / 7.0, 1.0e-9);
    EXPECT_NEAR(statistics.rms(), std::sqrt(29.0), 1.0e-9);
}

TEST(RunningStatistics, update_with_other_statistics_matches_merged_dataset) {
    math::RunningStatistics<double> first_half;
    for (const double sample : {2.0, 4.0, 4.0, 4.0}) {
        first_half.update(sample);
    }
    math::RunningStatistics<double> second_half;
    for (const double sample : {5.0, 5.0, 7.0, 9.0}) {
        second_half.update(sample);
    }
    first_half.update(second_half);
    EXPECT_EQ(first_half.num_samples(), 8u);
    EXPECT_DOUBLE_EQ(first_half.mean(), 5.0);
    EXPECT_DOUBLE_EQ(first_half.minimum(), 2.0);
    EXPECT_DOUBLE_EQ(first_half.maximum(), 9.0);
    EXPECT_NEAR(first_half.sum_of_square_differences(), 32.0, 1.0e-9);
    EXPECT_NEAR(first_half.population_variance(), 4.0, 1.0e-9);
}

TEST(RunningStatistics, update_with_self_throws) {
    math::RunningStatistics<double> statistics(5.0);
    EXPECT_THROW(statistics.update(statistics), std::runtime_error);
}

TEST(RunningStatistics, update_with_zero_sample_statistics_is_noop) {
    math::RunningStatistics<double> statistics(5.0);
    const math::RunningStatistics<double> empty;
    statistics.update(empty);
    EXPECT_EQ(statistics.num_samples(), 1u);
    EXPECT_DOUBLE_EQ(statistics.mean(), 5.0);
}

TEST(RunningStatistics, update_into_empty_statistics_copies_argument) {
    math::RunningStatistics<double> statistics;
    const math::RunningStatistics<double> other(5.0);
    statistics.update(other);
    EXPECT_EQ(statistics.num_samples(), 1u);
    EXPECT_DOUBLE_EQ(statistics.mean(), 5.0);
}

// Regression test: this reconstruction constructor used to compute sum_of_square_differences_ with (N-1) instead of
// N, and omit a term from sum_of_squares_, giving wrong population_variance()/rms() on the very next read. The
// constructed statistics here should exactly reproduce the values built up sample-by-sample above.
TEST(RunningStatistics, reconstruction_matches_incrementally_built_statistics) {
    const math::RunningStatistics<double> reconstructed(5.0, 4.0, 2.0, 9.0, 8);
    EXPECT_EQ(reconstructed.num_samples(), 8u);
    EXPECT_DOUBLE_EQ(reconstructed.mean(), 5.0);
    EXPECT_DOUBLE_EQ(reconstructed.minimum(), 2.0);
    EXPECT_DOUBLE_EQ(reconstructed.maximum(), 9.0);
    EXPECT_DOUBLE_EQ(reconstructed.sum(), 40.0);
    EXPECT_NEAR(reconstructed.sum_of_square_differences(), 32.0, 1.0e-9);
    EXPECT_NEAR(reconstructed.sum_of_squares(), 232.0, 1.0e-9);
    EXPECT_NEAR(reconstructed.rms(), std::sqrt(29.0), 1.0e-9);
    EXPECT_NEAR(reconstructed.sample_variance(), 32.0 / 7.0, 1.0e-9);
}

TEST(RunningStatistics, three_arg_reconstruction_sets_minimum_and_maximum_to_mean) {
    const math::RunningStatistics<double> reconstructed(5.0, 4.0, 8);
    EXPECT_DOUBLE_EQ(reconstructed.minimum(), 5.0);
    EXPECT_DOUBLE_EQ(reconstructed.maximum(), 5.0);
    EXPECT_NEAR(reconstructed.sum_of_square_differences(), 32.0, 1.0e-9);
}

TEST(RunningStatistics, zero_samples_with_nonzero_mean_throws) {
    EXPECT_THROW(math::RunningStatistics<double>(1.0, 0.0, 0), std::runtime_error);
}

TEST(RunningStatistics, zero_samples_with_zero_mean_and_variance_does_not_throw) {
    EXPECT_NO_THROW(math::RunningStatistics<double>(0.0, 0.0, 0));
}

TEST(RunningStatistics, one_sample_with_mismatched_minimum_throws) {
    EXPECT_THROW(math::RunningStatistics<double>(5.0, 0.0, 4.0, 5.0, 1), std::runtime_error);
}

TEST(RunningStatistics, one_sample_with_nonzero_variance_throws) {
    EXPECT_THROW(math::RunningStatistics<double>(5.0, 1.0, 1), std::runtime_error);
}

TEST(RunningStatistics, minimum_greater_than_maximum_throws) {
    EXPECT_THROW(math::RunningStatistics<double>(5.0, 0.0, 6.0, 4.0, 2), std::runtime_error);
}

// Regression test: the maximum-vs-mean invariant check used to compare maximum against the *sum* (num_samples *
// mean) rather than the mean itself, so it spuriously threw for entirely valid negative-mean data (e.g. samples
// {10, -5}: N=2, mean=2.5, sum=5, maximum=10 > sum, despite being valid).
TEST(RunningStatistics, valid_negative_mean_with_large_maximum_does_not_throw) {
    EXPECT_NO_THROW(math::RunningStatistics<double>(2.5, 0.0, -5.0, 10.0, 2));
}

TEST(RunningStatistics, maximum_less_than_mean_throws) {
    EXPECT_THROW(math::RunningStatistics<double>(0.0, 0.0, -2.0, -1.0, 2), std::runtime_error);
}

TEST(RunningStatisticsVector, default_construction_is_zero) {
    const math::RunningStatisticsVector<double> statistics(3);
    EXPECT_EQ(statistics.size(), 3u);
    EXPECT_TRUE(statistics.mean().isApprox(Eigen::Vector3d::Zero()));
    EXPECT_TRUE(statistics.population_variance().isApprox(Eigen::Vector3d::Zero()));
}

TEST(RunningStatisticsVector, update_with_eigen_vector) {
    math::RunningStatisticsVector<double> statistics(2);
    statistics.update(Eigen::Vector2d(1.0, 10.0));
    statistics.update(Eigen::Vector2d(3.0, 20.0));
    const Eigen::Vector2d mean = statistics.mean();
    EXPECT_DOUBLE_EQ(mean[0], 2.0);
    EXPECT_DOUBLE_EQ(mean[1], 15.0);
    const std::vector<std::size_t> num_samples = statistics.num_samples();
    EXPECT_EQ(num_samples[0], 2u);
    EXPECT_EQ(num_samples[1], 2u);
}

TEST(RunningStatisticsVector, update_with_std_vector) {
    math::RunningStatisticsVector<double> statistics(2);
    statistics.update(std::vector<double>{1.0, 10.0});
    statistics.update(std::vector<double>{3.0, 20.0});
    EXPECT_DOUBLE_EQ(statistics[0].mean(), 2.0);
    EXPECT_DOUBLE_EQ(statistics[1].mean(), 15.0);
}

TEST(RunningStatisticsVector, update_with_wrong_size_throws) {
    math::RunningStatisticsVector<double> statistics(2);
    EXPECT_THROW(statistics.update(std::vector<double>{1.0}), std::runtime_error);
    EXPECT_THROW(statistics.update(Eigen::Vector3d(1.0, 2.0, 3.0)), std::runtime_error);
}

TEST(RunningStatisticsVector, update_with_other_running_statistics_vector) {
    math::RunningStatisticsVector<double> statistics(2);
    statistics.update(std::vector<double>{1.0, 10.0});
    const std::vector<math::RunningStatistics<double>> other(2, math::RunningStatistics<double>(3.0));
    statistics.update(other);
    EXPECT_EQ(statistics[0].num_samples(), 2u);
    EXPECT_DOUBLE_EQ(statistics[0].mean(), 2.0);
    EXPECT_EQ(statistics[1].num_samples(), 2u);
    EXPECT_DOUBLE_EQ(statistics[1].mean(), 6.5);
}

TEST(RunningStatisticsVector, at_bounds_checked) {
    math::RunningStatisticsVector<double> statistics(2);
    EXPECT_NO_THROW(statistics.at(1));
    EXPECT_THROW(statistics.at(2), std::out_of_range);
}

TEST(RunningStatisticsVector, mutable_operator_bracket_updates_underlying_statistics) {
    math::RunningStatisticsVector<double> statistics(2);
    statistics[0].update(4.0);
    EXPECT_EQ(statistics.statistics()[0].num_samples(), 1u);
    EXPECT_DOUBLE_EQ(statistics.statistics()[0].mean(), 4.0);
}
