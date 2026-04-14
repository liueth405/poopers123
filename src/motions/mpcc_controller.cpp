#include "motions/mpcc_controller.hpp"
#include "acados/utils/types.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

// Parameter count verification
static_assert(
    MPCC_NP == 77,
    "Parameter count mismatch! Regenerate solver with generate_ocp.py");

MPCC_Controller::MPCC_Controller(double fk, double d_intake)
    : fk_(fk), d_intake_(d_intake), solve_time_ms_(0.0), sqp_iter_(0),
      kkt_norm_(0.0) {

  // Match Python dt_list: first 10 at 0.01s, remaining 15 at 0.02s (N=25)
  dt_list_.clear();
  for (int i = 0; i < 10 && i < MPCC_N; ++i)
    dt_list_.push_back(0.01);
  for (int i = 10; i < MPCC_N; ++i)
    dt_list_.push_back(0.02);

  capsule_ = mpcc_acados_create_capsule();

  int status =
      mpcc_acados_create_with_discretization(capsule_, MPCC_N, nullptr);
  if (status) {
    throw std::runtime_error(
        "mpcc_acados_create_with_discretization() failed with status " +
        std::to_string(status));
  }

  reset();
  setControlBounds(ControlBounds{});
}

MPCC_Controller::~MPCC_Controller() {
  if (capsule_) {
    mpcc_acados_free(capsule_);
    mpcc_acados_free_capsule(capsule_);
  }
}

void MPCC_Controller::reset() {
  ocp_nlp_config *nlp_config = mpcc_acados_get_nlp_config(capsule_);
  ocp_nlp_dims *nlp_dims = mpcc_acados_get_nlp_dims(capsule_);
  ocp_nlp_in *nlp_in = mpcc_acados_get_nlp_in(capsule_);
  ocp_nlp_out *nlp_out = mpcc_acados_get_nlp_out(capsule_);

  double x_init[MPCC_NX] = {0};
  double u_init[MPCC_NU] = {0};

  for (int i = 0; i <= MPCC_N; ++i) {
    ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, nlp_in, i, "x", x_init);
  }
  for (int i = 0; i < MPCC_N; ++i) {
    ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, nlp_in, i, "u", u_init);
  }
}

void MPCC_Controller::resetBuffer(const State &start_state) {
  ocp_nlp_config *nlp_config = mpcc_acados_get_nlp_config(capsule_);
  ocp_nlp_dims *nlp_dims = mpcc_acados_get_nlp_dims(capsule_);
  ocp_nlp_in *nlp_in = mpcc_acados_get_nlp_in(capsule_);
  ocp_nlp_out *nlp_out = mpcc_acados_get_nlp_out(capsule_);

  double x_init[MPCC_NX] = {start_state.x, start_state.y, start_state.theta,
                            0.0, start_state.s};
  double u_init[MPCC_NU] = {std::clamp(0.5, bounds_.v_min, bounds_.v_max), 0.0,
                            std::clamp(0.5, bounds_.vs_min, bounds_.vs_max)};

  for (int i = 0; i <= MPCC_N; ++i) {
    ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, nlp_in, i, "x", x_init);
  }
  for (int i = 0; i < MPCC_N; ++i) {
    ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, nlp_in, i, "u", u_init);
  }
}

void MPCC_Controller::setPath(const std::vector<PathVec2> &waypoints) {
  planner_.setWaypoints(waypoints);
}

void MPCC_Controller::setTerminalTargetVelocity(double v_end) {
  weights_.v_ref = v_end;
}

void MPCC_Controller::setIntakeOffset(double d_intake) { d_intake_ = d_intake; }

void MPCC_Controller::setControlBounds(const ControlBounds &b) {
  bounds_ = b;
  if (bounds_.v_min > bounds_.v_max || bounds_.w_min > bounds_.w_max ||
      bounds_.vs_min > bounds_.vs_max) {
    throw std::invalid_argument("MPCC_Controller::setControlBounds: min > max");
  }

  ocp_nlp_config *nlp_config = mpcc_acados_get_nlp_config(capsule_);
  ocp_nlp_dims *nlp_dims = mpcc_acados_get_nlp_dims(capsule_);
  ocp_nlp_in *nlp_in = mpcc_acados_get_nlp_in(capsule_);
  ocp_nlp_out *nlp_out = mpcc_acados_get_nlp_out(capsule_);

  double lbu[MPCC_NU] = {bounds_.v_min, bounds_.w_min, bounds_.vs_min};
  double ubu[MPCC_NU] = {bounds_.v_max, bounds_.w_max, bounds_.vs_max};

  for (int i = 0; i < MPCC_N; ++i) {
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i,
                                  "lbu", lbu);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i,
                                  "ubu", ubu);
  }
}

void MPCC_Controller::setWeights(const Weights &w) { weights_ = w; }

void MPCC_Controller::setSTrustRegion(double s_trust) {
  weights_.s_trust = s_trust;
}

void MPCC_Controller::updateParameters(const Control &prev_control,
                                       double current_s, double first_node_dt) {
  ocp_nlp_config *nlp_config = mpcc_acados_get_nlp_config(capsule_);
  ocp_nlp_dims *nlp_dims = mpcc_acados_get_nlp_dims(capsule_);
  ocp_nlp_out *nlp_out = mpcc_acados_get_nlp_out(capsule_);

  const double total_length = planner_.getTotalLength();
  static int last_ref_idx = 0;

  for (int i = 0; i < MPCC_N + 1; ++i) {
    double x_pred[5];
    ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, i, "x", x_pred);

    // Get predicted s for this node
    double s_node = std::clamp((double)x_pred[4], 0.0, total_length - 0.01);
    if (s_node < current_s)
      s_node = current_s;

    // Optimized parameter layout (77 total) aligned with generate_ocp.py:
    // p[0]: dt, p[1]: fk, p[2]: d_intake, p[3-5]: prev controls
    // p[6]: s_center (trust region), p[7]: q_pos, p[8]: q_heading
    // p[9]: q_cross, p[10]: q_vy, p[11]: q_progress, p[12]: w_dv
    // p[13]: w_dw, p[14]: w_dvs, p[15]: w_v, p[16]: v_ref
    // p[17]: q_pos_N, p[18]: q_heading_N, p[19]: q_cross_N, p[20]: q_vy_N
    // p[21]: kappa, p[22]: kappa_threshold, p[23]: s_trust
    // p[32-76]: Waypoints [x, y, vx, vy, s]

    double p[77] = {0};
    p[0] = (i < (int)dt_list_.size()) ? dt_list_[i] : dt_list_.back();
    if (i == 0 && first_node_dt > 0.0)
      p[0] = first_node_dt;

    p[1] = fk_;
    p[2] = d_intake_;
    p[3] = prev_control.v;
    p[4] = prev_control.w;
    p[5] = prev_control.vs;
    p[6] = s_node;

    p[7] = weights_.q_pos;
    p[8] = weights_.q_heading;
    p[9] = weights_.q_cross;
    p[10] = weights_.q_vy;
    p[11] = weights_.q_progress;
    p[12] = weights_.w_dv;
    p[13] = weights_.w_dw;
    p[14] = weights_.w_dvs;
    p[15] = weights_.w_v;
    p[16] = weights_.v_ref;

    p[17] = weights_.q_pos_N;
    p[18] = weights_.q_heading_N;
    p[19] = weights_.q_cross_N;
    p[20] = weights_.q_vy_N;

    // Estimate local curvature
    double eps_k = std::min(0.05, total_length - s_node);
    if (eps_k < 1e-4)
      eps_k = 1e-4;
    Pose pA = planner_.getPose(s_node, last_ref_idx);
    Pose pB =
        planner_.getPose(std::min(s_node + eps_k, total_length), last_ref_idx);
    double dth = pB.theta - pA.theta;
    while (dth > M_PI)
      dth -= 2.0 * M_PI;
    while (dth < -M_PI)
      dth += 2.0 * M_PI;
    double kappa_node = std::abs(dth) / eps_k;

    p[21] = kappa_node;
    p[22] = weights_.kappa_threshold;
    p[23] = weights_.s_trust;

    // WAYPOINT BUFFER (p[32..76]) - 9 waypoints for Hermite spline
    constexpr int NUM_WAYPOINTS = 9;
    double remaining = total_length - s_node;
    double lookahead = std::min(1.6, remaining);
    lookahead = std::max(lookahead, (NUM_WAYPOINTS - 1) * 0.02);

    for (int w_idx = 0; w_idx < NUM_WAYPOINTS; ++w_idx) {
      double sw = s_node + (w_idx * lookahead) / (NUM_WAYPOINTS - 1);
      sw = std::clamp(sw, 0.0, total_length);
      Pose pw = planner_.getPose(sw, last_ref_idx);

      int base = 32 + w_idx * 5;
      p[base + 0] = pw.x;
      p[base + 1] = pw.y;
      p[base + 2] = std::sin(pw.theta); // Tangent x (unit)
      p[base + 3] = std::cos(pw.theta); // Tangent y (unit)
      p[base + 4] = sw;
    }

    mpcc_acados_update_params(capsule_, i, p, 77);
  }
}

void MPCC_Controller::shiftPredictionHorizon() {
  ocp_nlp_config *nlp_config = mpcc_acados_get_nlp_config(capsule_);
  ocp_nlp_dims *nlp_dims = mpcc_acados_get_nlp_dims(capsule_);
  ocp_nlp_in *nlp_in = mpcc_acados_get_nlp_in(capsule_);
  ocp_nlp_out *nlp_out = mpcc_acados_get_nlp_out(capsule_);

  double x[MPCC_NX];
  double u[MPCC_NU];

  for (int i = 0; i < MPCC_N; ++i) {
    ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, i + 1, "x", x);
    ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, nlp_in, i, "x", x);
  }

  for (int i = 0; i < MPCC_N - 1; ++i) {
    ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, i + 1, "u", u);
    ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, nlp_in, i, "u", u);
  }
}

MPCC_Controller::Control MPCC_Controller::runStep(const State &current_state,
                                                  const Control &prev_control) {
  ocp_nlp_config *nlp_config = mpcc_acados_get_nlp_config(capsule_);
  ocp_nlp_dims *nlp_dims = mpcc_acados_get_nlp_dims(capsule_);
  ocp_nlp_in *nlp_in = mpcc_acados_get_nlp_in(capsule_);
  ocp_nlp_out *nlp_out = mpcc_acados_get_nlp_out(capsule_);
  ocp_nlp_solver *nlp_solver = mpcc_acados_get_nlp_solver(capsule_);

  double x0[MPCC_NX] = {current_state.x, current_state.y, current_state.theta,
                        current_state.vy, current_state.s};

  ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, 0, "lbx",
                                x0);
  ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, 0, "ubx",
                                x0);
  ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, nlp_in, 0, "x", x0);

  // Set control bounds
  double lbu[MPCC_NU] = {bounds_.v_min, bounds_.w_min, 0.0};
  double ubu[MPCC_NU] = {bounds_.v_max, bounds_.w_max, bounds_.vs_max};
  for (int i = 0; i < MPCC_N; ++i) {
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i,
                                  "lbu", lbu);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, nlp_out, i,
                                  "ubu", ubu);
  }

  updateParameters(prev_control, current_state.s, 0.0);

  int status = mpcc_acados_solve(capsule_);
  if (status != ACADOS_SUCCESS) {
    if (status != 2) {
      last_solve_.status = status;
      resetBuffer(current_state);
      return Control{0.0, 0.0, 0.0};
    }
  }

  double u_opt[MPCC_NU];
  ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, 0, "u", u_opt);

  ocp_nlp_get(nlp_solver, "time_tot", &solve_time_ms_);
  solve_time_ms_ *= 1000.0;
  ocp_nlp_get(nlp_solver, "sqp_iter", &sqp_iter_);
  ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, 0, "kkt_norm_inf", &kkt_norm_);

  last_solve_.status = ACADOS_SUCCESS;
  last_solve_.time_ms = solve_time_ms_;
  last_solve_.sqp_iter = sqp_iter_;
  last_solve_.kkt_norm = kkt_norm_;

  Control optimal_cmd;
  optimal_cmd.v = u_opt[0];
  optimal_cmd.w = u_opt[1];
  optimal_cmd.vs = u_opt[2];

  if (!std::isfinite(optimal_cmd.v) || !std::isfinite(optimal_cmd.w) ||
      !std::isfinite(optimal_cmd.vs)) {
    last_solve_.status = -999;
    resetBuffer(current_state);
    return Control{0.0, 0.0, 0.0};
  }

  shiftPredictionHorizon();
  return optimal_cmd;
}

double MPCC_Controller::projectSFromPosition(double x, double y,
                                             double s_hint) const {
  const double total = planner_.getTotalLength();
  if (total <= 1e-9)
    return 0.0;

  auto dist2_at = [&](double s) {
    Pose p = planner_.getPose(s);
    double dx = x - p.x;
    double dy = y - p.y;
    return dx * dx + dy * dy;
  };

  double window_back = std::min(0.3, s_hint * 0.15); // Larger backward window
  double window_fwd = (s_hint > 0.4) ? 0.25 : 0.40;  // Larger forward window for off-path recovery
  double s_min = std::clamp(s_hint - window_back, 0.0, total);
  double s_max = std::clamp(s_hint + window_fwd, 0.0, total);

  double best_s = s_hint;
  double best_d2 = dist2_at(std::clamp(best_s, 0.0, total));
  for (int i = 0; i <= 60; ++i) {
    double a = static_cast<double>(i) / 60.0;
    double s = s_min + a * (s_max - s_min);
    double d2 = dist2_at(s);
    if (d2 < best_d2) {
      best_d2 = d2;
      best_s = s;
    }
  }

  double r_min = std::clamp(best_s - 0.08, 0.0, total);
  double r_max = std::clamp(best_s + 0.08, 0.0, total);
  for (int i = 0; i <= 20; ++i) {
    double a = static_cast<double>(i) / 20.0;
    double s = r_min + a * (r_max - r_min);
    double d2 = dist2_at(s);
    if (d2 < best_d2) {
      best_d2 = d2;
      best_s = s;
    }
  }

  return std::clamp(best_s, 0.0, total);
}

double MPCC_Controller::getPathLength() const {
  return planner_.getTotalLength();
}
Pose MPCC_Controller::getPathPoseAtS(double s_m) const {
  return planner_.getPose(s_m);
}
double MPCC_Controller::getSolveTimeMs() const { return solve_time_ms_; }
int MPCC_Controller::getSQPIterations() const { return sqp_iter_; }
double MPCC_Controller::getKKTNorm() const { return kkt_norm_; }
