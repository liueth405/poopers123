#pragma once

#include "main.h"
#include <cstdint>
#include <cstdio>

/// Append-only CSV log of MPCC + estimate + chassis feedback on the V5 SD
/// card (same convention as sysid: paths like "/usd/...").
class MPCC_SdLog {
public:
  MPCC_SdLog() = default;
  MPCC_SdLog(const MPCC_SdLog &) = delete;
  MPCC_SdLog &operator=(const MPCC_SdLog &) = delete;

  /// Opens file. If path is nullptr or empty, uses "/usd/mpcc_<millis>.csv".
  /// Returns false if no SD or fopen fails.
  bool begin(const char *path_or_null);

  void end();

  bool is_open() const { return fp_ != nullptr; }

  void log_row(uint32_t t_ms, double path_len_m,
               // PF estimate (inches / deg except lateral in inches)
               double est_x_in, double est_y_in, double est_theta_deg,
               double est_lateral_in,
               // MPC state passed to solver (SI)
               double mpc_x_m, double mpc_y_m, double mpc_theta_rad,
               double mpc_vy_mps, double mpc_s_m,
               // Progress helpers (m)
               double s_proj_m, double s_err_clamped_m, double s_after_int_m,
               // MPC command (SI)
               double u_v_mps, double u_w_rads, double u_vs_mps,
               // Solver
               int solve_status, double solve_ms, int sqp_iter, double kkt_norm,
               // Path reference at mpc_s (SI)
               double ref_x_m, double ref_y_m, double ref_theta_rad,
               // Chassis EKF (before applying this step's velocity_control)
               double meas_v_ins, double meas_w_rads,
               // Commands sent to velocity_control this tick (SI)
               double cmd_v_mps, double cmd_w_rads);

private:
  FILE *fp_{nullptr};
  unsigned lines_{0};
};
