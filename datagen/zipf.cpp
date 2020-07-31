#include "zipf.h"

/**
 * Helper function that calculates log(1 + x) / x. A Taylor series expansion is
 * used if x is close to 0.
 * @param x A value larger than or equal to -1.
 * @return log(1 + x) / x.
 */
double helper1(double x) {
  if (std::abs(x) > 1e-8) {
    return std::log1p(x) / x;
  }
  return 1.0 - x * (0.5 - x * (1.0 / 3.0 - 0.25 * x));
}

/**
 * Helper function that calculates (exp(x) - 1) / x. A Taylor series expansion
 * is used if x is close to 0.
 * @param x Free parameter.
 * @return (exp(x) - 1) / x if x is nonzero, else 1.
 */
double helper2(double x) {
  if (std::abs(x) > 1e-8) {
    return std::expm1(x) / x;
  }
  return 1.0 + x * 0.5 * (1.0 + x / 3.0 * (1.0 + 0.25 * x));
}

ZipfDistribution::ZipfDistribution(int num_elements, double exponent)
    : num_elements_(num_elements > 0
                        ? num_elements
                        : throw std::logic_error(
                              "number of elements is not strictly positive: " +
                              std::to_string(num_elements))),
      exponent_(exponent > 0.0 ? exponent
                               : throw std::logic_error(
                                     "exponent is not strictly positive: " +
                                     std::to_string(exponent))),
      h_integral_x1_(h_integral(1.5) - 1),
      h_integral_num_elements_(h_integral(num_elements_ + 0.5)),
      s_(2 - h_integral_inv(h_integral(2.5) - h(2))) {}

double ZipfDistribution::h_integral(double x) const {
  double log_x = std::log(x);
  return helper2((1.0 - exponent_) * log_x) * log_x;
}

double ZipfDistribution::h(double x) const {
  return std::exp(-exponent_ * std::log(x));
}

double ZipfDistribution::h_integral_inv(double x) const {
  double t = x * (1.0 - exponent_);
  if (t < -1.0) {
    t = -1.0;
  }
  return std::exp(helper1(t) * x);
}
