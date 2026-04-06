#pragma once

#include "chassis.hpp"
#include "odom/odom.hpp"
#include "pros/rtos.hpp"
#include "utils.hpp"
#include <atomic>

/**
 * @brief Robot movement states.
 */
enum class MotionState { IDLE, DRIVING, TURNING, SETTLING };

/**
 * @brief Parameters for a Point-to-Point movement.
 */
struct PointTarget {
  double x, y;
  double max_v = 50.0, min_v = 2.0, slew = 1.0; // max_v in in/s, slew in in/s per tick
  bool reversed = false;
  PIDCoeffs linear_pid;
  double steering_gain = 1.5;
  double tolerance = 1.0; // inches
  double exit_turn_dist =
      4.0; // inches to target where we stop updating curvature
  double integral_limit = 0; // PID integral clamp
};

/**
 * @brief Parameters for a Turning movement.
 */
struct TurnTarget {
  double target_theta;
  enum class Direction { LEFT, RIGHT, NEAREST } direction = Direction::NEAREST;
  double max_v = 300.0, min_v = 5.0, slew = 6.0; // max_v in deg/s, slew in deg/s per tick
  PIDCoeffs angular_pid;
  double tolerance = 1.0; // degrees
  double integral_limit = 0;
};

/**
 * @brief MotionManager handles background motion control tasks and state
 * machine logic.
 */
class MotionManager {
public:
  /**
   * @brief Construct a new Motion Manager object.
   *
   * @param chassis Reference to initialized Chassis.
   * @param pf Reference to initialized ParticleFilter.
   */
  MotionManager(Chassis &chassis, ParticleFilter &pf);

  /**
   * @brief Start a Point-to-Point movement.
   */
  void moveTo(PointTarget target);

  /**
   * @brief Start a Turning movement to a specific angle.
   */
  void turnTo(TurnTarget target);

  /**
   * @brief Turn to face a specific coordinate.
   */
  void turnToPoint(double x, double y, TurnTarget::Direction dir,
                   PIDCoeffs angular_pid, double max_v = 300.0,
                   double slew = 6.0, double integral_limit = 0);

  /**
   * @brief Stop current motion and clear targets.
   */
  void cancelMotion();

  /**
   * @brief Check if target is reached.
   */
  bool isDone() const { return state == MotionState::IDLE; }

  /**
   * @brief Blocking wait functions.
   */
  void waitUntilPoint(double x, double y, double margin);
  void waitUntilAngle(double theta, double margin);
  void waitUntilPose(double x, double y, double theta, double pos_margin,
                     double theta_margin);
  void waitUntilDistance(double dist);

  /**
   * @brief Wait until linear speed matches condition
   * @param velocity Target speed in RPM (averaged across sides)
   * @param greaterThan True to wait for speed >= target, False for speed <=
   * target
   */
  void waitUntilLinearVelocity(double velocity, bool greaterThan = true);

  /**
   * @brief Wait until angular speed matches condition
   * @param omega Target angular speed in RPM (differential velocity)
   * @param greaterThan True to wait for speed >= target, False for speed <=
   * target
   */
  void waitUntilAngularVelocity(double omega, bool greaterThan = true);

  /**
   * @brief Unified high-frequency (10ms) update loop for Odom + Motions
   * This method is intended to run within a pros::Task
   */
  void loop();

private:
  Chassis &chassis;
  ParticleFilter &pf;

  std::atomic<MotionState> state{MotionState::IDLE};

  PointTarget current_point;
  TurnTarget current_turn;

  PID linear_controller = {{0, 0, 0}};
  PID angular_controller = {{0, 0, 0}};
  Slew v_slew = {500};
  Slew w_slew = {500};

  // Helper for angle wrapping
  double wrap_angle(double angle) const;
};
