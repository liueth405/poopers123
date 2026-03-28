#pragma once

#include <algorithm>
#include <cmath>

/**
 * @brief PID Coefficients.
 */
struct PIDCoeffs {
  double kP, kI, kD;
};

/**
 * @brief Simple Slew Rate helper to limit acceleration.
 */
class Slew {
public:
  double limit;
  double current = 0;

  /**
   * @brief Construct a new Slew object.
   *
   * @param limit Max change per call (e.g., max voltage change per 10ms).
   */
  Slew(double limit) : limit(limit) {}

  /**
   * @brief Update the slew-limited value.
   */
  double update(double target) {
    double delta = target - current;
    if (std::abs(delta) > limit) {
      current += (delta > 0 ? 1 : -1) * limit;
    } else {
      current = target;
    }
    return current;
  }

  /**
   * @brief Reset the slew value.
   */
  void reset(double val = 0) { current = val; }
};

/**
 * @brief Simple PID controller with Clamping Anti-integral Windup.
 */
class PID {
public:
  PIDCoeffs coeffs;
  double error = 0, prev_error = 0, integral = 0;
  double integral_limit = 0; // 0 means no limit

  PID(PIDCoeffs coeffs) : coeffs(coeffs) {}

  /**
   * @brief Update the PID controller.
   *
   * @param target Desired value.
   * @param current Current sensor value.
   * @return double Control output.
   */
  double update(double target, double current) {
    error = target - current;

    integral += error;
    // Anti-integral Windup (Clamping)
    if (integral_limit > 0) {
      if (integral > integral_limit)
        integral = integral_limit;
      if (integral < -integral_limit)
        integral = -integral_limit;
    }

    double derivative = error - prev_error;
    prev_error = error;
    return (error * coeffs.kP) + (integral * coeffs.kI) +
           (derivative * coeffs.kD);
  }

  /**
   * @brief Set the integral term limit.
   */
  void set_integral_limit(double limit) { integral_limit = limit; }

  /**
   * @brief Reset the PID state.
   */
  void reset() {
    error = 0;
    prev_error = 0;
    integral = 0;
  }
};
