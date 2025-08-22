#ifndef MATHBOX_IMPL_STATISTICS_HPP
#define MATHBOX_IMPL_STATISTICS_HPP

#include <cmath>
#include <cppbox/exceptions.hpp>
#include <stdexcept>

#include "mathbox/statistics.hpp"

namespace math {

template<std::floating_point Scalar_>
Statistics<Scalar_>::Statistics()
    : Statistics(static_cast<Scalar>(0), static_cast<Scalar>(0)) {}

template<std::floating_point Scalar_>
Statistics<Scalar_>::Statistics(const Scalar mean_, const Scalar population_variance_)
    : Statistics(mean_, population_variance_, mean_, mean_) {}

template<std::floating_point Scalar_>
Statistics<Scalar_>::Statistics(const Scalar mean_, const Scalar population_variance_, const Scalar minimum_,
        const Scalar maximum_)
    : mean_(mean_),
      population_variance_(population_variance_),
      minimum_(minimum_),
      maximum_(maximum_) {}

template<std::floating_point Scalar_>
inline auto Statistics<Scalar_>::maximum() const -> Scalar {
    return maximum_;
}

template<std::floating_point Scalar_>
inline auto Statistics<Scalar_>::maximum() -> Scalar& {
    return maximum_;
}

template<std::floating_point Scalar_>
inline auto Statistics<Scalar_>::mean() const -> Scalar {
    return mean_;
}

template<std::floating_point Scalar_>
inline auto Statistics<Scalar_>::mean() -> Scalar& {
    return mean_;
}

template<std::floating_point Scalar_>
inline auto Statistics<Scalar_>::minimum() const -> Scalar {
    return minimum_;
}

template<std::floating_point Scalar_>
inline auto Statistics<Scalar_>::minimum() -> Scalar& {
    return minimum_;
}

template<std::floating_point Scalar_>
inline auto Statistics<Scalar_>::population_stddev() const -> Scalar {
    return std::sqrt(population_variance());
}

template<std::floating_point Scalar_>
inline auto Statistics<Scalar_>::population_variance() const -> Scalar {
    return population_variance_;
}

template<std::floating_point Scalar_>
inline auto Statistics<Scalar_>::population_variance() -> Scalar& {
    return population_variance_;
}

template<std::floating_point Scalar_>
RunningStatistics<Scalar_>::RunningStatistics()
    : Statistics<Scalar_>(),
      num_samples_(0),
      sum_of_square_differences_(static_cast<Scalar>(0)),
      sum_(static_cast<Scalar_>(0)) {}

template<std::floating_point Scalar_>
RunningStatistics<Scalar_>::RunningStatistics(const Scalar sample)
    : Statistics<Scalar_>(sample, static_cast<Scalar>(0)),
      num_samples_(1),
      sum_of_square_differences_(static_cast<Scalar>(0)),
      sum_(sample) {}

template<std::floating_point Scalar_>
RunningStatistics<Scalar_>::RunningStatistics(const Scalar mean_, const Scalar population_variance_,
        const std::size_t num_samples_)
    : RunningStatistics<Scalar_>(mean_, population_variance_, mean_, mean_, num_samples_) {}

template<std::floating_point Scalar_>
RunningStatistics<Scalar_>::RunningStatistics(const Scalar mean_, const Scalar population_variance_,
        const Scalar minimum_, const Scalar maximum_, const std::size_t num_samples_)
    : Statistics<Scalar_>(mean_, population_variance_, minimum_, maximum_),
      num_samples_(num_samples_),
      sum_of_square_differences_(num_samples_ == 0 ? static_cast<Scalar>(0)
                                                   : static_cast<Scalar>(num_samples_ - 1) * population_variance_),
      sum_(static_cast<Scalar>(num_samples_) * mean_) {
    throw_if(num_samples_ == 0 && (mean_ != static_cast<Scalar>(0) || population_variance_ != static_cast<Scalar>(0)),
            "Cannot construct RunningStatistics with zero samples but non-zero mean or non-zero variance");
    if (num_samples_ == 1) [[unlikely]] {
        throw_if(mean_ != minimum_ || mean_ != maximum_,
                "Cannot construct RunningStatistics with one sample where minimum or maximum != mean");
        throw_if(population_variance_ != static_cast<Scalar>(0),
                "Cannot construct RunningStatistics with one sample but non-zero variance");
    }
    throw_if(minimum_ > maximum_, "Cannot construct RunningStatistics with minimum > maximum");
    throw_if(maximum_ > sum(),
            "Cannot construct RunningStatistics with maximum > sum (where sum = num_samples * mean)");
}

template<std::floating_point Scalar_>
inline std::size_t RunningStatistics<Scalar_>::num_samples() const {
    return num_samples_;
}

template<std::floating_point Scalar_>
inline auto RunningStatistics<Scalar_>::sample_stddev() const -> Scalar {
    return std::sqrt(sample_variance());
}

template<std::floating_point Scalar_>
inline auto RunningStatistics<Scalar_>::sample_variance() const -> Scalar {
    return num_samples() < 2 ? static_cast<Scalar>(0)
                             : sum_of_square_differences_ / static_cast<Scalar>(num_samples() - 1);
}

template<std::floating_point Scalar_>
inline auto RunningStatistics<Scalar_>::sum() const -> Scalar {
    return sum_;
}

template<std::floating_point Scalar_>
void RunningStatistics<Scalar_>::update(const Scalar sample) {
    ++num_samples_;
    sum_ += sample;
    if (num_samples_ == 1) {
        this->mean() = sample;
        this->minimum() = sample;
        this->maximum() = sample;
        return;
    }
    const Scalar previous_delta = sample - this->mean();
    this->mean() += previous_delta / static_cast<Scalar>(num_samples_);
    sum_of_square_differences_ += previous_delta * (sample - this->mean());
    this->population_variance() = sum_of_square_differences_ / num_samples_;
    this->minimum() = std::min(sample, this->minimum());
    this->maximum() = std::max(sample, this->maximum());
}

}

#if !MATHBOX_HEADER_ONLY
namespace math {

extern template class Statistics<double>;
extern template class RunningStatistics<double>;

}
#endif

#endif
