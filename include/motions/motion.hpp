#pragma once

#include "chassis.hpp"
#include "mpcc_controller.hpp"
#include "mpcc_sd_log.hpp"
#include "odom/odom.hpp"
#include "pathplanner.hpp"
#include "pros/rtos.hpp"
#include "utils.hpp"
#include <atomic>

/**
 * @brief Robot movement states.
 */
enum class MotionState {
  IDLE,
  DRIVING,
  TURNING,
  PRE_TURN_PATH,
  FOLLOWING_PATH
};

/**
 * @brief Parameters for a Point-to-Point movement.
 */
struct PointTarget {
  double x, y;
  double max_v = 50.0, min_v = 2.0,
         slew = 1.0; // max_v in in/s, slew in in/s per tick
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
  double max_v = 300.0, min_v = 5.0,
         slew = 6.0; // max_v in deg/s, slew in deg/s per tick
  PIDCoeffs angular_pid;
  double tolerance = 1.0; // degrees
  double integral_limit = 0;
};

/**
 * @brief Coordinate frame convention for path waypoints.
 */
enum class PathCoordFrame {
  FIELD_XY,   ///< Waypoints are (field_x, field_y) - standard Cartesian
  FWD_LATERAL ///< Waypoints are (forward, lateral) - robot-centric
};

/**
 * @brief Parameters for a Path-to-Path movement.
 */
struct PathTarget {
  std::vector<PathVec2> waypoints; // inches
  MPCC_Controller::Weights weights;
  /// MPC u bounds: v, vs in m/s; w in rad/s (same as acados/Python lbu/ubu).
  MPCC_Controller::ControlBounds mpc_control_bounds{};
  double tolerance;          // inches
  double v_target_end = 0.0; // Target terminal velocity (in/s)
  /// If true, append MPC debug CSV on SD (requires FAT32 card). See
  /// mpc_log_csv_path.
  bool log_mpc_to_sd = false;
  /// nullptr or "" -> "/usd/mpcc_<millis>.csv"; else exact path (e.g.
  /// "/usd/run.csv").
  const char *mpc_log_csv_path = nullptr;
  /// Coordinate frame of waypoints. If FWD_LATERAL, swaps X/Y to match field
  /// frame.
  PathCoordFrame coord_frame = PathCoordFrame::FIELD_XY;
  /// Distance from center of mass to intake, in inches.
  double intake_offset = 0.0;
  /// Slew rate limits (per 10ms tick) - prevents jittery oscillation on real
  /// hardware
  double max_dw_per_tick = 2.0; // rad/s per tick (default: 200 rad/s²)
  double max_dv_per_tick = 0.3; // m/s per tick (default: 30 m/s²)
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
  MotionManager(Chassis &chassis, ParticleFilter &pf,
                MPCC_Controller &mpcc_controller);

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
   * @brief Start a Path-to-Path movement.
   * @param target The target path to follow.
   */
  void followPath(PathTarget target);

  /**
   * @brief Wait until path is done.
   */
  void waitUntilPathDone();

  /**
   * @brief Unified high-frequency (10ms) update loop for Odom + Motions
   * This method is intended to run within a pros::Task
   */
  void loop();

private:
  Chassis &chassis;
  ParticleFilter &pf;
  MPCC_Controller &mpcc_controller;

  std::atomic<MotionState> state{MotionState::IDLE};

  PointTarget current_point;
  TurnTarget current_turn;

  PID linear_controller = {{0, 0, 0}};
  PID angular_controller = {{0, 0, 0}};
  Slew v_slew = {500};
  Slew w_slew = {500};

  PathTarget current_path;
  double current_mpcc_s = 0.0;
  MPCC_Controller::Control prev_mpcc_u = {0.0, 0.0, 0.0};
  MPCC_SdLog mpc_sd_log_;

  // Helper for angle wrapping
  double wrap_angle(double angle) const;
};
