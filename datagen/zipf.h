#ifndef JOINSREVISITED_DATAGEN_ZIPF_H_
#define JOINSREVISITED_DATAGEN_ZIPF_H_

#include <cmath>
#include <random>

/**
 * @brief Implementation of the <a
 * href="https://en.wikipedia.org/wiki/Zipf's_law">Zipf distribution</a>.
 * @details This is a direct port of Apache Common's <a
 * href="https://github.com/apache/commons-rng/blob/6a1b0c16090912e8fc5de2c1fb5bd8490ac14699/commons-rng-sampling/src/main/java/org/apache/commons/rng/sampling/distribution/RejectionInversionZipfSampler.java">RejectionInverseZipfSampler</a>.
 * It uses rejection inversion sampling for a discrete, bounded Zipf
 * distribution that is based on the method described by Wolfgang Hörmann and
 * Gerhard Derflinger in <a
 * href="https://dl.acm.org/citation.cfm?id=235029"><i>Rejection-inversion to
 * generate variates from monotone discrete distributions</i></a>, ACM
 * Transactions on Modeling and Computer Simulation (TOMACS) 6.3 (1996):
 * 169-184.
 */
class ZipfDistribution {
public:
  /**
   * @param num_elements Number of elements.
   * @param exponent Exponent.
   * @throws std::logic_error if num_elements and exponent are not strictly
   * positive.
   */
  ZipfDistribution(int num_elements, double exponent);

  /**
   * Generate a random number from the Zipf distribution.
   * @tparam R A random number generator type (e.g., std::mt19937).
   * @param rng A random number generator.
   * @return A random number from the Zipf distribution.
   */
  template <class R> int operator()(R &rng) {
    // The paper describes an algorithm for exponents larger than 1 (Algorithm
    // ZRI). The original method uses
    //   H(x) = (v + x)^(1 - q) / (1 - q)
    // as the integral of the hat function. This function is undefined for q =
    // 1, which is the reason for the limitation of the exponent. If instead the
    // integral function
    //   H(x) = ((v + x)^(1 - q) - 1) / (1 - q)
    // is used, for which a meaningful limit exists for q = 1, the method works
    // for all positive exponents. The following implementation uses v = 0 and
    // generates integral number in the range [1, num_elements]. This is
    // different to the original method where v is defined to be positive and
    // numbers are taken from [0, i_max]. This explains why the implementation
    // looks slightly different.

    while (true) {
      double u =
          h_integral_num_elements_ +
          real_distribution_(rng) * (h_integral_x1_ - h_integral_num_elements_);
      // u is uniformly distributed in (h_integral_x1_,
      // h_integral_num_elements_].

      double x = h_integral_inv(u);
      int k = (int)std::round(x);

      // Limit k to the range [1, num_elements_] if it would be outside due
      // to numerical inaccuracies.
      if (k < 1) {
        k = 1;
      } else if (k > num_elements_) {
        k = num_elements_;
      }

      // Here, the distribution of k is given by:
      //   P(k = 1) = C * (h_integral(1.5) - h_integral_x1_) = C
      //   P(k = m) = C * (h_integral(m + 1/2) - h_integral(m - 1/2)) for m >= 2
      // where C = 1 / (h_integral_num_elements_ - h_integral_x1_)

      if (k - x <= s_ || u >= h_integral(k + 0.5) - h(k)) {
        // Case k = 1:
        //
        //   The right inequality is always true, because replacing k by 1 gives
        //   u >= h_integral(1.5) - h(1) = h_integral_x1_ and u is taken from
        //   (h_integral_x1_, h_integral_num_elements_].
        //
        //   Therefore, the acceptance rate for k = 1 is P(accepted | k = 1) = 1
        //   and the probability that 1 is returned as random value is
        //   P(k = 1 and accepted) = P(accepted | k = 1) * P(k = 1) = C = C /
        //   1^exponent
        //
        // Case k >= 2:
        //
        //   The left inequality (k - x <= s) is just a short cut
        //   to avoid the more expensive evaluation of the right inequality
        //   (u >= h_integral(k + 0.5) - h(k)) in many cases.
        //
        //   If the left inequality is true, the right inequality is also true:
        //     Theorem 2 in the paper is valid for all positive exponents,
        //     because the requirements h'(x) = -exponent/x^(exponent + 1) < 0
        //     and
        //     (-1/h_integral_inverse'(x))'' = (1+1/exponent) * x^(1/exponent-1)
        //     >= 0 are both fulfilled. Therefore, f(x) = x -
        //     h_integral_inv(h_integral(x + 0.5) - h(x)) is a non-decreasing
        //     function. If k - x <= s holds, k - x <= s + f(k) - f(2) is
        //     obviously also true which is equivalent to -x <=
        //     -h_integral_inv(h_integral(k + 0.5) - h(k)), -h_integral_inv(u)
        //     <= -h_integral_inv(h_integral(k + 0.5) - h(k)), and finally u >=
        //     h_integral(k + 0.5) - h(k).
        //
        //   Hence, the right inequality determines the acceptance rate:
        //   P(accepted | k = m) = h(m) / (hIntegrated(m+1/2) -
        //   hIntegrated(m-1/2)) The probability that m is returned is given by
        //   P(k = m and accepted) = P(accepted | k = m) * P(k = m) = C * h(m) =
        //   C / m^exponent.
        //
        // In both cases the probabilities are proportional to the probability
        // mass function of the Zipf distribution.

        return k;
      }
    }
  }

private:
  /// The number of elements. This is the upper bound (inclusive) on the random
  /// number that can be generated.
  int num_elements_;

  /// The exponent parameter of the Zipf distribution.
  double exponent_;

  /// h_integral(1.5) - 1.
  double h_integral_x1_;

  /// h_integral(num_elements_ + 0.5).
  double h_integral_num_elements_;

  /// 2 - h_integral_inv(h_integral(2.5) - h(2)).
  double s_;

  /// A real distribution between 0.0 and 1.0.
  std::uniform_real_distribution<double> real_distribution_;

  /**
   * H(x), the integral function of h(x).
   * @param x Free parameter.
   * @return H(x).
   */
  [[nodiscard]] double h_integral(double x) const;

  /**
   * h(x) = 1 / x^exponent_.
   * @param x Free parameter.
   * @return h(x).
   */
  [[nodiscard]] double h(double x) const;

  /**
   * The inverse function of H(x).
   * @param x Free parameter.
   * @return y for which H(y) = x.
   */
  [[nodiscard]] double h_integral_inv(double x) const;
};

#endif // JOINSREVISITED_DATAGEN_ZIPF_H_
