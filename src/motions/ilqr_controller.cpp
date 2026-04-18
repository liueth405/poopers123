#include "motions/ilqr_controller.hpp"
#include "pros/rtos.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>

// Indexing helpers
#define IDX_K(k, i, j)                                                         \
  ((k) * ILQR_Controller::NU * ILQR_Controller::NX +                           \
   (i) * ILQR_Controller::NX + (j))
#define IDX_K_U(k, i) ((k) * ILQR_Controller::NU + (i))

ILQR_Controller::ILQR_Controller(double fk, double d_intake)
    : fk_(fk), d_intake_(d_intake) {
  std::memset(x_, 0, sizeof(x_));
  std::memset(y_, 0, sizeof(y_));
  std::memset(th_, 0, sizeof(th_));
  std::memset(vy_, 0, sizeof(vy_));
  std::memset(s_, 0, sizeof(s_));
  std::memset(v_, 0, sizeof(v_));
  std::memset(w_, 0, sizeof(w_));
  std::memset(vs_, 0, sizeof(vs_));
  std::memset(K_, 0, sizeof(K_));
  std::memset(k_, 0, sizeof(k_));
}

void ILQR_Controller::setPath(const std::vector<PathWaypoint> &waypoints) {
  planner_.setWaypoints(waypoints);
}

void ILQR_Controller::saturateControl(Control &u) {
  // Clamp vs only (v and w are normalized separately)
  u.vs = std::clamp(u.vs, bounds_.vs_min, bounds_.vs_max);
}

// Normalize v and w to respect differential drive wheel speed limits
// For diff drive: v_left = v - w*track/2, v_right = v + w*track/2
// Scale (v, w) proportionally if either wheel exceeds limits
void ILQR_Controller::normalizeWheelSpeeds(Control &u) {
  double half_track = bounds_.track_width * 0.5;
  double v_left = u.v - u.w * half_track;
  double v_right = u.v + u.w * half_track;

  // Find max violation (check both min and max)
  double max_violation = 1.0;
  if (v_left > bounds_.max_wheel_speed)
    max_violation = std::max(max_violation, v_left / bounds_.max_wheel_speed);
  if (v_right > bounds_.max_wheel_speed)
    max_violation = std::max(max_violation, v_right / bounds_.max_wheel_speed);
  if (v_left < bounds_.min_wheel_speed)
    max_violation = std::max(max_violation, bounds_.min_wheel_speed / v_left);
  if (v_right < bounds_.min_wheel_speed)
    max_violation = std::max(max_violation, bounds_.min_wheel_speed / v_right);

  if (max_violation > 1.0) {
    double scale = 1.0 / max_violation;
    u.v *= scale;
    u.w *= scale;
  }
}

// Drift model from generate_ocp.py:
// f_k_dt = fk * dt
// centripetal = v * w * dt
// cent_clamped = f_k_dt * tanh(centripetal / (f_k_dt + eps))
// vy_sign = tanh(100.0 * vy) // Smooth sign
// vy_next = vy - cent_clamped - vy_sign * f_k_dt
// Differentiable friction-drift model matching odom.hpp simplified physics
// v_pre = vy - (v * w * dt)
// v_next = v_pre - f_k_dt * tanh(C * v_pre)
double ILQR_Controller::computeDriftModel(double vy, double v, double w) const {
  double f_k_dt = fk_ * dt_;
  double v_pre = vy - (v * w * dt_);
  // Use a sharp tanh to simulate resistive friction without changing sign
  return v_pre - f_k_dt * std::tanh(10.0 * v_pre);
}

void ILQR_Controller::computeDynamics(int k) {
  // Compass heading: theta=0 is North, CW positive
  // dx/dt = v*sin(th) + vy*cos(th)
  // dy/dt = v*cos(th) - vy*sin(th)

  double th_mid = th_[k] + 0.5 * w_[k] * dt_;
  double sin_mid = std::sin(th_mid);
  double cos_mid = std::cos(th_mid);

  x_[k + 1] = x_[k] + (v_[k] * sin_mid + vy_[k] * cos_mid) * dt_;
  y_[k + 1] = y_[k] + (v_[k] * cos_mid - vy_[k] * sin_mid) * dt_;
  th_[k + 1] = th_[k] + w_[k] * dt_;

  // Drift model from generate_ocp.py
  vy_[k + 1] = computeDriftModel(vy_[k], v_[k], w_[k]);

  s_[k + 1] = s_[k] + vs_[k] * dt_;

  // Update trust region center
  s_center_[k + 1] = s_[k + 1];
}

void ILQR_Controller::computeCostDerivatives(int k, bool is_terminal,
                                             double *lx, double *lu,
                                             double lxx[NX][NX],
                                             double luu[NU][NU],
                                             double lux[NU][NX]) {
  double s_ref = std::clamp(s_[k], 0.0, planner_.getTotalLength());
  Pose ref = planner_.getPose(s_ref);

  // Position error
  double dx = x_[k] - ref.x;
  double dy = y_[k] - ref.y;
  double pos_error = std::sqrt(dx * dx + dy * dy + 1e-6);

  // Path tangent and normal (compass: sin=dx/ds, cos=dy/ds)
  double path_tan_x = std::sin(ref.theta);
  double path_tan_y = std::cos(ref.theta);
  double path_norm_x = std::cos(ref.theta); // Left of forward
  double path_norm_y = -std::sin(ref.theta);

  // Cross-track error (signed, perpendicular to path)
  double cross_track = dx * path_norm_x + dy * path_norm_y;

  // Heading error (continuous formulation)
  double costh = std::cos(th_[k]);
  double sinth = std::sin(th_[k]);
  double cosref = std::cos(ref.theta);
  double sinref = std::sin(ref.theta);

  // Gradient: 2*sin(th - ref)
  double heading_grad = 2.0 * (sinth * cosref - costh * sinref);
  // Hessian: 2*cos(th - ref)
  double heading_hess = 2.0 * (costh * cosref + sinth * sinref);

  // Select weights based on terminal vs stage
  double w_pos = is_terminal ? weights_.q_pos_N : weights_.q_pos;
  double w_heading = is_terminal ? weights_.q_heading_N : weights_.q_heading;
  double w_cross = is_terminal ? weights_.q_cross_N : weights_.q_cross;
  double w_vy = is_terminal ? weights_.q_vy_N : weights_.q_vy;

  // ========== STATE GRADIENT lx ==========
  // Position error gradient (linear cost: w_pos * pos_error)
  lx[0] = w_pos * dx / pos_error;
  lx[1] = w_pos * dy / pos_error;

  // Cross-track error gradient: d/dx (0.5 * w_cross * cross^2) = w_cross *
  // cross * d_cross/dx
  lx[0] += w_cross * cross_track * path_norm_x;
  lx[1] += w_cross * cross_track * path_norm_y;

  // Heading gradient
  lx[2] = w_heading * heading_grad;

  // Lateral velocity gradient
  lx[3] = w_vy * vy_[k];

  // Progress gradient (stage only)
  lx[4] = is_terminal ? 0.0 : -weights_.q_progress;
  // ========== CONTROL GRADIENT lu (stage only) ==========
  if (!is_terminal) {
    double dv_prev = v_[k] - v_prev_;
    double dw_prev = w_[k] - w_prev_;
    double dvs_prev = vs_[k] - vs_prev_;

    // Regularization: w_v * v^2 + w_w * w^2 + w_vs * vs^2
    lu[0] = 2.0 * weights_.w_v * v_[k] + weights_.w_dv * dv_prev -
            weights_.q_progress;
    lu[1] = 2.0 * weights_.w_w * w_[k] + weights_.w_dw * dw_prev;
    lu[2] = 2.0 * weights_.w_vs * vs_[k] + weights_.w_dvs * dvs_prev -
            weights_.q_progress;

    // Optional: add tracking cost if v_ref > 0
    if (fabs(weights_.v_ref) > 0.01) {
      lu[0] += 2.0 * weights_.w_v * (v_[k] - weights_.v_ref);
      lu[2] += 2.0 * weights_.w_vs * (vs_[k] - weights_.v_ref);
    }
  }

  // ========== STATE HESSIAN lxx ==========
  std::memset(lxx, 0, NX * NX * sizeof(double));

  // Position Hessian (from cross-track)
  lxx[0][0] = w_cross * path_norm_x * path_norm_x;
  lxx[1][1] = w_cross * path_norm_y * path_norm_y;
  lxx[0][1] = lxx[1][0] = w_cross * path_norm_x * path_norm_y;

  // Heading Hessian
  lxx[2][2] = w_heading * heading_hess;

  // Lateral velocity Hessian
  lxx[3][3] = w_vy;

  // ========== CONTROL HESSIAN luu (stage only) ==========
  if (!is_terminal) {
    std::memset(luu, 0, NU * NU * sizeof(double));
    luu[0][0] = weights_.w_v + weights_.w_dv;
    luu[1][1] = weights_.w_w + weights_.w_dw;
    luu[2][2] = weights_.w_vs + weights_.w_dvs;

    // Optional: add tracking Hessian if v_ref > 0
    if (weights_.v_ref > 0.01) {
      luu[0][0] += weights_.w_v;
      luu[2][2] += weights_.w_vs;
    }
  }

  // ========== CROSS HESSIAN lux ==========
  std::memset(lux, 0, NU * NX * sizeof(double));

  // ========== TRUST REGION (stage only) ==========
  if (!is_terminal) {
    double s_upper = s_center_[k] + weights_.s_trust;
    double s_lower = s_center_[k] - weights_.s_trust;

    double s_viol = 0.0;
    if (s_[k] > s_upper) {
      s_viol = s_[k] - s_upper;
    } else if (s_[k] < s_lower) {
      s_viol = s_lower - s_[k];
    }

    if (s_viol > 0.0) {
      // Add trust region cost: penalty * s_viol^2
      // Gradient: 2 * penalty * s_viol * sign(violation)
      double s_grad = 2.0 * weights_.s_trust_penalty * s_viol;
      if (s_[k] > s_upper) {
        lx[4] += s_grad; // Positive gradient when over upper bound
      } else {
        lx[4] -= s_grad; // Negative gradient when under lower bound
      }
      // Hessian contribution
      lxx[4][4] += 2.0 * weights_.s_trust_penalty;
    }
  }
}

void ILQR_Controller::backwardPass() {
  double lx[NX], lu[NU];
  double lxx[NX][NX], luu[NU][NU], lux[NU][NX];

  // Terminal cost
  computeCostDerivatives(N, true, lx, lu, lxx, luu, lux);

  double Vx[NX];
  double Vxx[NX][NX];

  for (int i = 0; i < NX; i++) {
    Vx[i] = lx[i];
    for (int j = 0; j < NX; j++)
      Vxx[i][j] = lxx[i][j];
  }

  // Regularize terminal Hessian
  for (int i = 0; i < NX; i++)
    Vxx[i][i] += mu_;

  // Backward recursion
  for (int k = N - 1; k >= 0; k--) {
    computeCostDerivatives(k, false, lx, lu, lxx, luu, lux);

    // Dynamics Jacobians
    double th_mid = th_[k] + 0.5 * w_[k] * dt_;
    double sin_mid = std::sin(th_mid);
    double cos_mid = std::cos(th_mid);

    // State Jacobian A (5x5)
    double A[NX][NX] = {};
    A[0][0] = 1.0;
    A[1][1] = 1.0;
    A[2][2] = 1.0;
    A[4][4] = 1.0;

    A[0][2] = (v_[k] * cos_mid - vy_[k] * sin_mid) * dt_;
    A[1][2] = -(v_[k] * sin_mid + vy_[k] * cos_mid) * dt_;
    A[0][3] = cos_mid * dt_;
    A[1][3] = -sin_mid * dt_;

    // Drift model Jacobian A[3][3] = d(vy_next)/dvy
    double f_k_dt = fk_ * dt_;
    double v_pre = vy_[k] - (v_[k] * w_[k] * dt_);
    double th_val = std::tanh(10.0 * v_pre);
    double dv_pre_dvy = 1.0;
    double df_dvpre = 1.0 - f_k_dt * 10.0 * (1.0 - th_val * th_val);
    A[3][3] = df_dvpre * dv_pre_dvy;

    // Control Jacobian B (5x3)
    double B[NX][NU] = {};
    B[0][0] = sin_mid * dt_;
    B[1][0] = cos_mid * dt_;
    B[2][1] = dt_;
    B[4][2] = dt_;

    // dvy_next / dv and dvy_next / dw
    double dv_pre_dv = -w_[k] * dt_;
    double dv_pre_dw = -v_[k] * dt_;
    B[3][0] = df_dvpre * dv_pre_dv;
    B[3][1] = df_dvpre * dv_pre_dw;

    // Compute Q-function
    double Qx[NX] = {};
    for (int i = 0; i < NX; i++) {
      Qx[i] = lx[i];
      for (int j = 0; j < NX; j++)
        Qx[i] += A[j][i] * Vx[j];
    }

    double Qu[NU] = {};
    for (int i = 0; i < NU; i++) {
      Qu[i] = lu[i];
      for (int j = 0; j < NX; j++)
        Qu[i] += B[j][i] * Vx[j];
    }

    double Vxx_A[NX][NX] = {};
    for (int i = 0; i < NX; i++) {
      for (int j = 0; j < NX; j++) {
        for (int p = 0; p < NX; p++)
          Vxx_A[i][j] += Vxx[i][p] * A[p][j];
      }
    }

    double Qxx[NX][NX] = {};
    for (int i = 0; i < NX; i++) {
      for (int j = 0; j < NX; j++) {
        Qxx[i][j] = lxx[i][j];
        for (int p = 0; p < NX; p++)
          Qxx[i][j] += A[p][i] * Vxx_A[p][j];
      }
    }

    double Vxx_B[NX][NU] = {};
    for (int i = 0; i < NX; i++) {
      for (int j = 0; j < NU; j++) {
        for (int p = 0; p < NX; p++)
          Vxx_B[i][j] += Vxx[i][p] * B[p][j];
      }
    }

    double Quu[NU][NU] = {};
    for (int i = 0; i < NU; i++) {
      for (int j = 0; j < NU; j++) {
        Quu[i][j] = luu[i][j];
        for (int p = 0; p < NX; p++)
          Quu[i][j] += B[p][i] * Vxx_B[p][j];
      }
    }

    double Qux[NU][NX] = {};
    for (int i = 0; i < NU; i++) {
      for (int j = 0; j < NX; j++) {
        Qux[i][j] = lux[i][j];
        for (int p = 0; p < NX; p++)
          Qux[i][j] += B[p][i] * Vxx_A[p][j];
      }
    }

    // Regularize Quu
    for (int i = 0; i < NU; i++)
      Quu[i][i] += mu_;

    // Direct 3x3 inverse
    double det = Quu[0][0] * (Quu[1][1] * Quu[2][2] - Quu[1][2] * Quu[2][1]) -
                 Quu[0][1] * (Quu[1][0] * Quu[2][2] - Quu[1][2] * Quu[2][0]) +
                 Quu[0][2] * (Quu[1][0] * Quu[2][1] - Quu[1][1] * Quu[2][0]);

    if (std::abs(det) < 1e-12) {
      mu_ = std::min(mu_ * 10.0, mu_max_);
      continue;
    }

    double inv_det = 1.0 / det;
    double Quu_inv[NU][NU];
    Quu_inv[0][0] = (Quu[1][1] * Quu[2][2] - Quu[1][2] * Quu[2][1]) * inv_det;
    Quu_inv[0][1] = (Quu[0][2] * Quu[2][1] - Quu[0][1] * Quu[2][2]) * inv_det;
    Quu_inv[0][2] = (Quu[0][1] * Quu[1][2] - Quu[0][2] * Quu[1][1]) * inv_det;
    Quu_inv[1][0] = (Quu[1][2] * Quu[2][0] - Quu[1][0] * Quu[2][2]) * inv_det;
    Quu_inv[1][1] = (Quu[0][0] * Quu[2][2] - Quu[0][2] * Quu[2][0]) * inv_det;
    Quu_inv[1][2] = (Quu[0][2] * Quu[1][0] - Quu[0][0] * Quu[1][2]) * inv_det;
    Quu_inv[2][0] = (Quu[1][0] * Quu[2][1] - Quu[1][1] * Quu[2][0]) * inv_det;
    Quu_inv[2][1] = (Quu[0][1] * Quu[2][0] - Quu[0][0] * Quu[2][1]) * inv_det;
    Quu_inv[2][2] = (Quu[0][0] * Quu[1][1] - Quu[0][1] * Quu[1][0]) * inv_det;

    // K = -Quu_inv * Qux
    for (int i = 0; i < NU; i++) {
      for (int j = 0; j < NX; j++) {
        K_[IDX_K(k, i, j)] = 0.0;
        for (int p = 0; p < NU; p++)
          K_[IDX_K(k, i, j)] -= Quu_inv[i][p] * Qux[p][j];
      }
    }

    // k = -Quu_inv * Qu
    for (int i = 0; i < NU; i++) {
      k_[IDX_K_U(k, i)] = 0.0;
      for (int p = 0; p < NU; p++)
        k_[IDX_K_U(k, i)] -= Quu_inv[i][p] * Qu[p];
    }

    // Expected improvement
    double dV = 0.0;
    for (int i = 0; i < NU; i++) {
      for (int j = 0; j < NU; j++)
        dV -= 0.5 * Qu[i] * Quu_inv[i][j] * Qu[j];
    }

    if (dV < 0)
      mu_ = std::min(mu_ * 10.0, mu_max_);
    else
      mu_ = std::max(mu_ / 10.0, mu_min_);

    // Update value function
    for (int i = 0; i < NX; i++) {
      Vx[i] = Qx[i];
      for (int j = 0; j < NU; j++)
        Vx[i] -= K_[IDX_K(k, j, i)] * Qu[j];
    }

    double KtQuu[NX][NU] = {};
    for (int i = 0; i < NX; i++) {
      for (int j = 0; j < NU; j++) {
        for (int p = 0; p < NU; p++)
          KtQuu[i][j] += K_[IDX_K(k, p, i)] * Quu[p][j];
      }
    }

    for (int i = 0; i < NX; i++) {
      for (int j = 0; j < NX; j++) {
        Vxx[i][j] = Qxx[i][j];
        for (int p = 0; p < NU; p++)
          Vxx[i][j] += KtQuu[i][p] * K_[IDX_K(k, p, j)];
      }
    }
  }
}

void ILQR_Controller::forwardPass(const State &x0, double alpha) {
  x_[0] = x0.x;
  y_[0] = x0.y;
  th_[0] = x0.theta;
  vy_[0] = x0.vy;
  s_[0] = x0.s;
  s_center_[0] = x0.s;

  if (alpha < 1.0) {
    std::memcpy(x_nom_, x_, sizeof(x_));
    std::memcpy(y_nom_, y_, sizeof(y_));
    std::memcpy(th_nom_, th_, sizeof(th_));
    std::memcpy(vy_nom_, vy_, sizeof(vy_));
    std::memcpy(s_nom_, s_, sizeof(s_));
  }

  for (int k = 0; k < N; k++) {
    double dx[NX] = {};
    if (alpha < 1.0) {
      dx[0] = x_[k] - x_nom_[k];
      dx[1] = y_[k] - y_nom_[k];
      dx[2] = th_[k] - th_nom_[k];
      dx[3] = vy_[k] - vy_nom_[k];
      dx[4] = s_[k] - s_nom_[k];
    }

    double dv_fb = 0.0, dw_fb = 0.0, dvs_fb = 0.0;
    for (int j = 0; j < NX; j++) {
      dv_fb += K_[IDX_K(k, 0, j)] * dx[j];
      dw_fb += K_[IDX_K(k, 1, j)] * dx[j];
      dvs_fb += K_[IDX_K(k, 2, j)] * dx[j];
    }

    double v_new = v_[k] + alpha * (k_[IDX_K_U(k, 0)] + dv_fb);
    double w_new = w_[k] + alpha * (k_[IDX_K_U(k, 1)] + dw_fb);
    double vs_new = vs_[k] + alpha * (k_[IDX_K_U(k, 2)] + dvs_fb);

    v_[k] = v_new;
    w_[k] = w_new;
    vs_[k] = std::clamp(vs_new, bounds_.vs_min, bounds_.vs_max);
    // Normalize to respect wheel speed limits
    Control temp{v_[k], w_[k], vs_[k]};
    normalizeWheelSpeeds(temp);
    v_[k] = temp.v;
    w_[k] = temp.w;

    computeDynamics(k);
  }
}

double ILQR_Controller::computeCost() {
  double cost = 0.0;
  double total_length = planner_.getTotalLength();

  for (int k = 0; k <= N; k++) {
    double s_ref = std::clamp(s_[k], 0.0, total_length);
    Pose ref = planner_.getPose(s_ref);

    double dx = x_[k] - ref.x;
    double dy = y_[k] - ref.y;
    double pos_error = std::sqrt(dx * dx + dy * dy + 1e-6);

    double heading_cost = (std::cos(th_[k]) - std::cos(ref.theta)) *
                              (std::cos(th_[k]) - std::cos(ref.theta)) +
                          (std::sin(th_[k]) - std::sin(ref.theta)) *
                              (std::sin(th_[k]) - std::sin(ref.theta));

    // Cross-track error
    double path_norm_x = std::cos(ref.theta);
    double path_norm_y = -std::sin(ref.theta);
    double cross_track = dx * path_norm_x + dy * path_norm_y;

    bool is_terminal = (k == N);
    double w_pos = is_terminal ? weights_.q_pos_N : weights_.q_pos;
    double w_heading = is_terminal ? weights_.q_heading_N : weights_.q_heading;
    double w_cross = is_terminal ? weights_.q_cross_N : weights_.q_cross;
    double w_vy = is_terminal ? weights_.q_vy_N : weights_.q_vy;
    // Stage cost
    cost += w_pos * pos_error;
    cost += 0.5 * w_heading * heading_cost;
    cost += 0.5 * w_cross * cross_track * cross_track;
    cost += 0.5 * w_vy * vy_[k] * vy_[k];

    if (k < N) {
      // Control costs (regularization + smoothness)
      cost += 0.5 * weights_.w_v * v_[k] * v_[k];
      cost += 0.5 * weights_.w_w * w_[k] * w_[k];
      cost += 0.5 * weights_.w_vs * vs_[k] * vs_[k];
      cost += 0.5 * weights_.w_dv * (v_[k] - v_prev_) * (v_[k] - v_prev_);
      cost += 0.5 * weights_.w_dw * (w_[k] - w_prev_) * (w_[k] - w_prev_);
      cost += 0.5 * weights_.w_dvs * (vs_[k] - vs_prev_) * (vs_[k] - vs_prev_);

      // Progress reward
      cost -= weights_.q_progress * vs_[k];

      // Trust region
      double s_upper = s_center_[k] + weights_.s_trust;
      double s_lower = s_center_[k] - weights_.s_trust;
      double s_viol = 0.0;
      if (s_[k] > s_upper)
        s_viol = s_[k] - s_upper;
      else if (s_[k] < s_lower)
        s_viol = s_lower - s_[k];
      cost += weights_.s_trust_penalty * s_viol * s_viol;
    }
  }

  return cost;
}

ILQR_Controller::Control ILQR_Controller::runStep(const State &current_state,
                                                  const Control &prev_control) {
  uint32_t start_time = pros::millis();

  v_prev_ = prev_control.v;
  w_prev_ = prev_control.w;
  vs_prev_ = prev_control.vs;

  x_[0] = current_state.x;
  y_[0] = current_state.y;
  th_[0] = current_state.theta;
  vy_[0] = current_state.vy;
  s_[0] = current_state.s;
  s_center_[0] = current_state.s;

  // Initial rollout
  for (int k = 0; k < N; k++) {
    double v_init = (weights_.v_ref > 0.01) ? weights_.v_ref : 0.5;
    v_[k] = v_init;
    w_[k] = 0.0;
    vs_[k] = std::clamp(v_init, bounds_.vs_min, bounds_.vs_max);
    Control temp{v_[k], w_[k], vs_[k]};
    normalizeWheelSpeeds(temp);
    v_[k] = temp.v;
    w_[k] = temp.w;
    computeDynamics(k);
  }

  double cost_old = computeCost();

  // Main iLQR loop
  for (int iter = 0; iter < 10; iter++) {
    iterations_used_ = iter + 1;
    backwardPass();

    double alpha = 1.0;
    double cost_new = cost_old;
    bool improved = false;

    double x_saved[N + 1], y_saved[N + 1], th_saved[N + 1];
    double vy_saved[N + 1], s_saved[N + 1];
    double v_saved[N], w_saved[N], vs_saved[N];

    std::memcpy(x_saved, x_, sizeof(x_));
    std::memcpy(y_saved, y_, sizeof(y_));
    std::memcpy(th_saved, th_, sizeof(th_));
    std::memcpy(vy_saved, vy_, sizeof(vy_));
    std::memcpy(s_saved, s_, sizeof(s_));
    std::memcpy(v_saved, v_, sizeof(v_));
    std::memcpy(w_saved, w_, sizeof(w_));
    std::memcpy(vs_saved, vs_, sizeof(vs_));

    for (int ls = 0; ls < 5; ls++) {
      forwardPass(current_state, alpha);
      cost_new = computeCost();

      if (cost_new < cost_old - 1e-4 * std::abs(cost_old)) {
        improved = true;
        break;
      }

      std::memcpy(x_, x_saved, sizeof(x_));
      std::memcpy(y_, y_saved, sizeof(y_));
      std::memcpy(th_, th_saved, sizeof(th_));
      std::memcpy(vy_, vy_saved, sizeof(vy_));
      std::memcpy(s_, s_saved, sizeof(s_));
      std::memcpy(v_, v_saved, sizeof(v_));
      std::memcpy(w_, w_saved, sizeof(w_));
      std::memcpy(vs_, vs_saved, sizeof(vs_));

      alpha *= 0.5;
    }

    if (improved) {
      cost_old = cost_new;
      mu_ = std::max(mu_ / 10.0, mu_min_);
    } else {
      mu_ = std::min(mu_ * 10.0, mu_max_);
    }

    if (std::abs(cost_old - cost_new) < 1e-4 * std::abs(cost_old) && improved)
      break;
  }

  solve_time_ms_ = static_cast<double>(pros::millis() - start_time);

  Control u;
  u.v = v_[0];
  u.w = w_[0];
  u.vs = vs_[0];
  normalizeWheelSpeeds(u); // Scale v,w to respect wheel speed limits
  return u;
}

double ILQR_Controller::projectSFromPosition(double x, double y,
                                             double s_hint) const {
  // Find closest point on path to (x, y)
  // Start from s_hint and search nearby
  double total_length = planner_.getTotalLength();
  if (total_length < 1e-6)
    return 0.0;

  // Search in a window around s_hint
  double best_s = s_hint;
  double best_dist_sq = 1e20;

  double search_radius = 0.5; // 0.5m search window
  double s_start = std::max(0.0, s_hint - search_radius);
  double s_end = std::min(total_length, s_hint + search_radius);

  // Coarse search
  int coarse_steps = 20;
  double coarse_step = (s_end - s_start) / coarse_steps;
  for (int i = 0; i <= coarse_steps; ++i) {
    double s = s_start + i * coarse_step;
    Pose p = planner_.getPose(s);
    double dx = x - p.x;
    double dy = y - p.y;
    double dist_sq = dx * dx + dy * dy;
    if (dist_sq < best_dist_sq) {
      best_dist_sq = dist_sq;
      best_s = s;
    }
  }

  // Fine search around best coarse point
  double fine_radius = coarse_step * 2;
  double fine_start = std::max(0.0, best_s - fine_radius);
  double fine_end = std::min(total_length, best_s + fine_radius);
  int fine_steps = 20;
  double fine_step = (fine_end - fine_start) / fine_steps;

  for (int i = 0; i <= fine_steps; ++i) {
    double s = fine_start + i * fine_step;
    Pose p = planner_.getPose(s);
    double dx = x - p.x;
    double dy = y - p.y;
    double dist_sq = dx * dx + dy * dy;
    if (dist_sq < best_dist_sq) {
      best_dist_sq = dist_sq;
      best_s = s;
    }
  }

  return best_s;
}
