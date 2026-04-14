#include "motions/mpcc_sd_log.hpp"
#include "pros/misc.hpp"
#include <cstdio>

namespace {
constexpr unsigned kFlushEvery = 25;
}

bool MPCC_SdLog::begin(const char *path_or_null) {
  end();
  if (pros::usd::is_installed() != 1) {
    std::printf("[MPCC_SdLog] SD card not installed\n");
    return false;
  }

  char path[64];
  if (path_or_null != nullptr && path_or_null[0] != '\0') {
    std::snprintf(path, sizeof(path), "%s", path_or_null);
  } else {
    std::snprintf(path, sizeof(path), "/usd/mpcc_%lu.csv",
                  static_cast<unsigned long>(pros::millis()));
  }

  fp_ = std::fopen(path, "w");
  if (!fp_) {
    std::printf("[MPCC_SdLog] fopen failed: %s\n", path);
    return false;
  }

  std::fprintf(fp_,
               "t_ms,path_len_m,"
               "est_x_in,est_y_in,est_theta_deg,est_lat_in,"
               "mpc_x_m,mpc_y_m,mpc_th_rad,mpc_vy_mps,mpc_s_m,"
               "s_proj_m,s_err_clamped_m,s_after_int_m,"
               "u_v_mps,u_w_rads,u_vs_mps,"
               "solve_status,solve_ms,sqp_iter,kkt,"
               "ref_x_m,ref_y_m,ref_th_rad,"
               "meas_v_ins,meas_w_rads,cmd_v_mps,cmd_w_rads\n");
  std::fflush(fp_);
  lines_ = 0;
  return true;
}

void MPCC_SdLog::end() {
  if (fp_) {
    std::fflush(fp_);
    std::fclose(fp_);
    fp_ = nullptr;
  }
  lines_ = 0;
}

void MPCC_SdLog::log_row(
    uint32_t t_ms, double path_len_m, double est_x_in, double est_y_in,
    double est_theta_deg, double est_lateral_in, double mpc_x_m,
    double mpc_y_m, double mpc_theta_rad, double mpc_vy_mps, double mpc_s_m,
    double s_proj_m, double s_err_clamped_m, double s_after_int_m,
    double u_v_mps, double u_w_rads, double u_vs_mps, int solve_status,
    double solve_ms, int sqp_iter, double kkt_norm, double ref_x_m,
    double ref_y_m, double ref_theta_rad, double meas_v_ins, double meas_w_rads,
    double cmd_v_mps, double cmd_w_rads) {
  if (!fp_)
    return;

  std::fprintf(
      fp_,
      "%u,%.6f,"
      "%.4f,%.4f,%.4f,%.4f,"
      "%.6f,%.6f,%.6f,%.6f,%.6f,"
      "%.6f,%.6f,%.6f,"
      "%.6f,%.6f,%.6f,"
      "%d,%.4f,%d,%.6f,"
      "%.6f,%.6f,%.6f,"
      "%.6f,%.6f,%.6f,%.6f\n",
      static_cast<unsigned>(t_ms), path_len_m, est_x_in, est_y_in,
      est_theta_deg, est_lateral_in, mpc_x_m, mpc_y_m, mpc_theta_rad,
      mpc_vy_mps, mpc_s_m, s_proj_m, s_err_clamped_m, s_after_int_m, u_v_mps,
      u_w_rads, u_vs_mps, solve_status, solve_ms, sqp_iter, kkt_norm, ref_x_m,
      ref_y_m, ref_theta_rad, meas_v_ins, meas_w_rads, cmd_v_mps, cmd_w_rads);

  ++lines_;
  if (lines_ % kFlushEvery == 0)
    std::fflush(fp_);
}
