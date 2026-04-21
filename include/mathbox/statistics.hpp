#ifndef MATHBOX_STATISTICS_HPP
#define MATHBOX_STATISTICS_HPP

#include <Eigen/Core>
#include <concepts>
#include <cstdint>
#include <vector>

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

    /**
     * @brief Get number of samples.
     *
     * @return Scalar
     */
    std::size_t num_samples() const;

    /**
     * @brief Root mean square:
     * \f[
     *      \text{RMS} = \sqrt{\frac{1}{N} \sum_i^N{x_i^2}}
     * \f]
     *
     * @return Scalar
     */
    Scalar rms() const;

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
     * @brief Get the sum of squares.
     *
     * @return Scalar
     */
    Scalar sum_of_squares() const;

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

private:
    std::size_t num_samples_;
    Scalar sum_of_square_differences_;
    Scalar sum_;
    Scalar sum_of_squares_;
};

template<std::floating_point Scalar_>
class RunningStatisticsVector {
public:
    using Scalar = Scalar_;

    explicit RunningStatisticsVector(const std::size_t size);

    const std::vector<RunningStatistics<Scalar>>& statistics() const;

    std::vector<RunningStatistics<Scalar>>& statistics();

    const RunningStatistics<Scalar>& at(const std::size_t i) const;

    RunningStatistics<Scalar>& at(const std::size_t i);

    const RunningStatistics<Scalar>& operator[](const std::size_t i) const;

    RunningStatistics<Scalar>& operator[](const std::size_t i);

    Eigen::Vector<Scalar, Eigen::Dynamic> rms() const;

    std::size_t size() const;

    void update(const Eigen::Vector<Scalar, Eigen::Dynamic>& samples);

    void update(const std::vector<Scalar>& samples);

    void update(const std::vector<RunningStatistics<Scalar>>& statistics__);

private:
    std::vector<RunningStatistics<Scalar>> statistics_;
};

}

#include "mathbox/impl/statistics.hpp"

#endif
