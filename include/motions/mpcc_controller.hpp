#pragma once

#include "motions/pathplanner.hpp"

#include <array>
#include <vector>

// Acados
#include "acados_c/ocp_nlp_interface.h"
#include "mpcc/acados_solver_mpcc.h"

class MPCC_Controller {
public:
  struct State {
    double x;
    double y;
    double theta;
    double vy;
    double s;
  };

  struct Control {
    double v;
    double w;
    double vs;
  };

  /// Box constraints on u = [v, w, vs]
  /// v, vs in m/s; w in rad/s.
  struct ControlBounds {
    double v_min = 0.0;
    double v_max = 1.5;
    double w_min = -2.5;
    double w_max = 2.5;
    double vs_min = 0.0;
    double vs_max = 2.0;
  };

  // Runtime-configurable weights (simplified vs MPCC)
  // Uses Euclidean distance + heading error instead of contouring/lag decomposition
  struct Weights {
    double q_pos = 10.0;      // Position error weight (Euclidean distance)
    double q_heading = 5.0;   // Heading error weight
    double q_vy = 4.0;        // Lateral velocity penalty
    double q_progress = 2.0;  // Progress reward (linear, encourages forward motion)
        double q_cross = 8.0;      // Cross-track error weight (sideways from path)

    double w_v = 6.0;         // Velocity tracking weight
    double v_ref = 0.8;       // Target velocity [m/s]

    // Control rate penalties
    double w_dv = 5.0;
    double w_dw = 5.0;
    double w_dvs = 3.0;

    // Terminal weights
    double q_pos_N = 20.0;
    double q_heading_N = 10.0;
    double q_cross_N = 16.0;   // Terminal cross-track weight
        double q_vy_N = 8.0;

    // s trust region
    double s_trust = 0.3;     // Half-width [m]

    // Curvature-based adaptation
    double kappa_threshold = 4.0;
  };

  MPCC_Controller(double fk = 0.3, double d_intake = 0.25);
  ~MPCC_Controller();

  void reset();
  void resetBuffer(const State &start_state);

  void setPath(const std::vector<PathVec2> &waypoints);

  void setTerminalTargetVelocity(double v_end);
  void setIntakeOffset(double d_intake);
  double getIntakeOffset() const { return d_intake_; }
  void setControlBounds(const ControlBounds &bounds);
  void setWeights(const Weights &weights);
  void setSTrustRegion(double s_trust);
  void setFrictionCoefficient(double fk) { fk_ = fk; }

  Control runStep(const State &current_state, const Control &prev_control);
  double projectSFromPosition(double x, double y, double s_hint) const;
  double getPathLength() const;
  Pose getPathPoseAtS(double s_m) const;

  double getSolveTimeMs() const;
  int getSQPIterations() const;
  double getKKTNorm() const;

  struct SolveInfo {
    int status = -1;
    double time_ms = 0.0;
    int sqp_iter = 0;
    double kkt_norm = 0.0;
  };
  SolveInfo getLastSolveInfo() const { return last_solve_; }

  void updateParameters(const Control &prev_control, double current_s,
                        double first_node_dt = 0.0);

private:
  mpcc_solver_capsule *capsule_;
  PathPlanner planner_;

  std::vector<double> dt_list_;
  double fk_;
  double d_intake_;
  Weights weights_;

  ControlBounds bounds_;

  double solve_time_ms_;
  int sqp_iter_;
  double kkt_norm_;
  SolveInfo last_solve_;

  void shiftPredictionHorizon();
};
