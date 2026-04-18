#pragma once

#include "motions/pathplanner.hpp"
#include <cstdint>
#include <cstring>

/**
 * Hand-written iLQR controller for path following.
 * Matches generate_ocp.py implementation.
 *
 * States: [x, y, theta, vy, s] (5D) - compass heading
 * Controls: [v, w, vs] (3D)
 *
 * Compass heading: theta=0 is North (positive Y), CW positive
 * - sin(theta) = dx/ds (x is East)
 * - cos(theta) = dy/ds (y is North)
 *
 * Features:
 * - Cross-track error (perpendicular to path)
 * - Drift model: vy_next = vy - fk*dt*tanh(v*w*dt/(fk*dt)) - sign(vy)*fk*dt
 * - Trust region soft constraint on s deviation
 * - Separate terminal weights
 * - Wheel speed based control constraints (differential drive)
 */

class ILQR_Controller {
public:
  // State: x, y, theta, lateral velocity vy, path progress s
  struct State {
    double x;
    double y;
    double theta;
    double vy;
    double s;
  };

  // Control: forward velocity, angular velocity, path progress rate
  struct Control {
    double v;
    double w;
    double vs;
  };

  // Weights for cost function (separate stage vs terminal)
  struct Weights {
    // Stage weights
    double q_pos = 10.0;     // Position error (linear Euclidean distance)
    double q_heading = 5.0;  // Heading error weight
    double q_cross = 8.0;    // Cross-track error weight
    double q_vy = 4.0;       // Lateral velocity penalty
    double q_progress = 2.0; // Progress reward (incentive to move forward)

    // Control regularization (penalize magnitude, not tracking)
    double w_v = 1.0;  // Velocity magnitude penalty (v^2)
    double w_w = 0.5;  // Angular velocity magnitude penalty (w^2)
    double w_vs = 0.5; // Path progress rate magnitude penalty (vs^2)

    // Control smoothness (penalize changes)
    double w_dv = 5.0;  // Acceleration penalty
    double w_dw = 5.0;  // Angular acceleration penalty
    double w_dvs = 3.0; // Progress rate smoothness penalty

    // Target velocity (optional - set > 0 to enable tracking)
    double v_ref = 0.0; // If > 0, adds tracking cost to this velocity

    // Terminal weights (higher for stability)
    double q_pos_N = 20.0;
    double q_heading_N = 10.0;
    double q_cross_N = 16.0;
    double q_vy_N = 8.0;

    // Trust region
    double s_trust = 0.15;          // Trust region half-width for progress
    double s_trust_penalty = 100.0; // Trust region violation penalty
  };

  // Control bounds (wheel speed based for differential drive)
  struct Bounds {
    double min_wheel_speed = -2.0; // m/s min individual wheel speed (reverse)
    double max_wheel_speed = 2.0;  // m/s max individual wheel speed
    double track_width = 0.35;     // meters between wheels
    double vs_min = 0.0;
    double vs_max = 2.0;
  };

  // Horizon length
  static constexpr int N = 50; // 50 steps @ 20ms = 1000ms horizon

  // State and control dimensions
  static constexpr int NX = 5;
  static constexpr int NU = 3;

  ILQR_Controller(double fk = 0.3, double d_intake = 0.25);
  ~ILQR_Controller() = default;

  // Set path
  void setPath(const std::vector<PathWaypoint> &waypoints);

  // Set weights and bounds
  void setWeights(const Weights &w) { weights_ = w; }
  void setBounds(const Bounds &b) { bounds_ = b; }
  void setTerminalVelocity(double v_ref) { weights_.v_ref = v_ref; }
  void setIntakeOffset(double offset) { d_intake_ = offset; }
  void setFrictionCoeff(double fk) { fk_ = fk; }
  void setTrustRegion(double width, double penalty) {
    weights_.s_trust = width;
    weights_.s_trust_penalty = penalty;
  }

  // Run one iteration
  Control runStep(const State &current_state, const Control &prev_control);

  // Get path info
  double getPathLength() const { return planner_.getTotalLength(); }
  Pose getPathPoseAtS(double s) const { return planner_.getPose(s); }

  // Stats
  double projectSFromPosition(double x, double y, double s_hint) const;
  double getSolveTimeMs() const { return solve_time_ms_; }
  int getIterations() const { return iterations_used_; }

private:
  PathPlanner planner_;
  Weights weights_;
  Bounds bounds_;
  double d_intake_ = 0.25; // Distance from center to intake point
  double fk_ = 0.3;        // Friction coefficient for drift model

  // Trajectory storage (stack-allocated for speed)
  double x_[N + 1];
  double y_[N + 1];
  double th_[N + 1];
  double vy_[N + 1];
  double s_[N + 1];
  double v_[N];
  double w_[N];
  double vs_[N];

  // Nominal trajectory for deviation calculation
  double x_nom_[N + 1];
  double y_nom_[N + 1];
  double th_nom_[N + 1];
  double vy_nom_[N + 1];
  double s_nom_[N + 1];

  // Previous control for smoothness
  double v_prev_ = 0.0;
  double w_prev_ = 0.0;
  double vs_prev_ = 0.0;

  // Feedback gains K (N x NU x NX) and feedforward k (N x NU)
  double K_[N * NU * NX];
  double k_[N * NU];

  // Trust region center (updated each solve)
  double s_center_[N + 1];

  // Tuning
  double dt_ = 0.02; // 20ms timestep
  double mu_ = 1e-4; // Levenberg-Marquardt regularization
  double mu_min_ = 1e-6;
  double mu_max_ = 1e-2;

  // Stats
  double solve_time_ms_ = 0.0;
  int iterations_used_ = 0;

  // Internal methods
  void forwardPass(const State &x0, double alpha);
  void backwardPass();
  double computeCost();
  void computeDynamics(int k);
  void computeCostDerivatives(int k, bool is_terminal, double *lx, double *lu,
                              double lxx[NX][NX], double luu[NU][NU],
                              double lux[NU][NX]);
  void saturateControl(Control &u); // Clamp vs only
  void
  normalizeWheelSpeeds(Control &u); // Scale v,w to respect wheel speed limits

  // Drift model helper
  double computeDriftModel(double vy, double v, double w) const;
};
