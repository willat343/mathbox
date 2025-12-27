#ifndef MATHBOX_STATISTICS_HPP
#define MATHBOX_STATISTICS_HPP

#include <concepts>
#include <cstdint>

namespace math {

template<std::floating_point Scalar_>
class Statistics {
public:
    /**
     * @brief Floating-point scalar type.
     *
     */
    using Scalar = Scalar_;

    /**
     * @brief Construct zero-initialised statistics.
     *
     */
    explicit Statistics();

    /**
     * @brief Construct statistics with known mean and variance.
     *
     * @param mean_
     * @param population_variance_
     */
    explicit Statistics(const Scalar mean_, const Scalar population_variance_);

    /**
     * @brief Construct statistics with known mean and variance, minimum and maximum.
     *
     * @param mean_
     * @param population_variance_
     * @param minimum_
     * @param maximum_
     */
    explicit Statistics(const Scalar mean_, const Scalar population_variance_, const Scalar minimum_,
            const Scalar maximum_);

    /**
     * @brief Get maximum.
     *
     * @return Scalar
     */
    Scalar maximum() const;

    /**
     * @brief Get mutable maximum.
     *
     * @return Scalar&
     */
    Scalar& maximum();

    /**
     * @brief Get mean.
     *
     * @return Scalar
     */
    Scalar mean() const;

    /**
     * @brief Get mutable mean.
     *
     * @return Scalar&
     */
    Scalar& mean();

    /**
     * @brief Get minimum.
     *
     * @return Scalar
     */
    Scalar minimum() const;

    /**
     * @brief Get mutable minimum.
     *
     * @return Scalar&
     */
    Scalar& minimum();

    /**
     * @brief Compute population stddev.
     *
     * @return Scalar
     */
    Scalar population_stddev() const;

    /**
     * @brief Get population variance
     *
     * @return Scalar
     */
    Scalar population_variance() const;

    /**
     * @brief Get mutable population variance
     *
     * @return Scalar&
     */
    Scalar& population_variance();

private:
    Scalar mean_;
    Scalar population_variance_;
    Scalar minimum_;
    Scalar maximum_;
};

template<std::floating_point Scalar_>
class RunningStatistics : public Statistics<Scalar_> {
public:
    /**
     * @brief Floating-point scalar type
     *
     */
    using Scalar = Scalar_;

    /**
     * @brief Construct running statistics with no samples (zero-mean zero-variance).
     *
     */
    explicit RunningStatistics();

    /**
     * @brief Construct running statistics from first sample.
     *
     * @param sample
     */
    explicit RunningStatistics(const Scalar sample);

    /**
     * @brief Construct statistics with known (current) mean and variance from a certain number of samples.
     *
     * @param num_samples_
     * @param mean_
     * @param population_variance_
     */
    explicit RunningStatistics(const Scalar mean_, const Scalar population_variance_, const std::size_t num_samples_);

    /**
     * @brief Construct statistics with known (current) mean and variance from a certain number of samples, and known
     * minimum, and maximum.
     *
     * @param mean_
     * @param population_variance_
     * @param minimum_
     * @param maximum_
     * @param num_samples_
     */
    explicit RunningStatistics(const Scalar mean_, const Scalar population_variance_, const Scalar minimum_,
            const Scalar maximum_, const std::size_t num_samples_);

    using Statistics<Scalar>::mean;

    /**
     * @brief Get number of samples.
     *
     * @return Scalar
     */
    std::size_t num_samples() const;

    using Statistics<Scalar>::population_variance;

    /**
     * @brief Compute sample standard deviation.
     *
     * @return Scalar
     */
    Scalar sample_stddev() const;

    /**
     * @brief Compute sample variance.
     *
     * @return Scalar
     */
    Scalar sample_variance() const;

    /**
     * @brief Get the sum of samples. A running sum is maintained, more numerically stable than num_samples() * mean().
     *
     * @return Scalar
     */
    Scalar sum() const;

    /**
     * @brief Get the sum of square differences.
     *
     * @return Scalar
     */
    Scalar sum_of_square_differences() const;

    /**
     * @brief Update statistics with a new sample. For the mean and variances, Welford's Algorithm is used.
     *
     * @param sample
     */
    void update(const Scalar sample);

    /**
     * @brief Update statistics with another statistics objects (representing samples from the same population), also
     * known as a Chan–Golub–LeVeque merge.
     *
     * @param statistics
     */
    void update(const RunningStatistics& statistics);

protected:
    /**
     * @brief Remove mutable variance access.
     *
     * @return Scalar&
     */
    inline Scalar& mean() {
        return Statistics<Scalar>::mean();
    }

    /**
     * @brief Remove mutable population variance access.
     *
     * @return Scalar&
     */
    inline Scalar& population_variance() {
        return Statistics<Scalar>::population_variance();
    }

private:
    std::size_t num_samples_;
    Scalar sum_of_square_differences_;
    Scalar sum_;
};

}

#include "mathbox/impl/statistics.hpp"

#endif
