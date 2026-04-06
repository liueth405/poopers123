#pragma once

#include "api.h"
#include "blasfeo/blasfeo.h"

class SystemIdentification;
#include <cmath>
#include <memory>
#include <vector>

/**
 * @brief Simple 2x2 Matrix for dynamics calculations.
 * Can be used as a fallback or for easy assignment.
 */
struct Mat2 {
  double data[4] = {0, 0, 0, 0};
  double &operator()(int r, int c) { return data[r * 2 + c]; }
  const double &operator()(int r, int c) const { return data[r * 2 + c]; }
  void setZero() {
    for (int i = 0; i < 4; ++i)
      data[i] = 0;
  }
  Mat2 &operator<<(double v) {
    data[idx++] = v;
    if (idx >= 4)
      idx = 0;
    return *this;
  }
  Mat2 &operator,(double v) {
    data[idx++] = v;
    if (idx >= 4)
      idx = 0;
    return *this;
  }

private:
  int idx = 0;
};

struct Vec2 {
  double data[2] = {0, 0};
  double &operator()(int i) { return data[i]; }
  const double &operator()(int i) const { return data[i]; }
  static Vec2 Zero() { return {0, 0}; }
};

/**
 * @brief Parameters defining the motor's physical properties.
 */
struct MotorParams {
  double R;    // resistance [ohm]
  double Kt;   // torque constant [N*m/A]
  double Ke;   // back-EMF constant [V*s/rad] (≈ Kt in SI)
  double Kv;   // velocity constant [rad/s/V] (≈ 1/Ke)
  double Imax; // current limit [A]
  double Vmax; // voltage limit [V]

  void normalize() {
    if (Ke <= 0 && Kv > 0)
      Ke = 1.0 / Kv;
    else if (Kv <= 0 && Ke > 0)
      Kv = 1.0 / Ke;
  }
};

/**
 * @brief Gains for the Forward/Gains dynamics model.
 */
struct GainsFG {
  double alpha_v;
  double alpha_w;
  double lambda_v;
  double lambda_w;
  double max_accel_wall_v;
  double max_accel_wall_w;
  double kS_v;
  double kS_w;
};

/**
 * @brief Parameters for the Extended Kalman Filter (EKF) velocity estimator.
 */
struct EKFConfig {
  double Q_v = 0.05; // Process noise variance for linear velocity
  double Q_w = 0.05; // Process noise variance for angular velocity
  double R_v_enc =
      0.02; // Measurement noise variance for encoder linear velocity
  double R_w_enc =
      0.02; // Measurement noise variance for encoder angular velocity
  double R_w_imu = 0.005; // Measurement noise variance for IMU angular velocity
  double slip_tolerance =
      9.0; // Mahalanobis distance threshold to reject slip (Chi-square)
};

/**
 * @brief Parameters for the Super-Twisting Sliding Mode Controller (STSMC).
 */
struct STSMCParams {
  double Lambda[2] = {1.2, 1.2};
  double K1[2] = {6.5, 6.5};
  double K2[2] = {12.0, 12.0};
  double s_eps = 0.01;
  double z_max[2] = {12.0, 12.0};
  double Imax_per_side = 5.0;
};

/**
 * @brief Chassis class to handle robot movement and sensor data.
 */
class Chassis {
  friend class SystemIdentification;
public:
  struct MotorData {
    double position = 0;
    double velocity = 0; // Estimated RPM
    double est_vel = 0;  // Filtered velocity state (combined or side-based)
    bool connected = false;
    double voltage = 0;
  };

  /**
   * @brief Construct a new Chassis object.
   */
  Chassis(std::vector<int> left_ports, std::vector<int> right_ports,
          std::vector<int> imu_ports, std::vector<int> distance_ports,
          double gear_ratio, double wheelDiameter, double trackWidth,
          pros::v5::MotorGears cartridge, EKFConfig config = EKFConfig());

  /**
   * @brief Destroy the Chassis object, freeing blasfeo memory.
   */
  ~Chassis();

  void tank(int left_voltage, int right_voltage);

  double get_rotation();
  double get_left_position();
  double get_right_position();
  double get_left_velocity();
  double get_right_velocity();
  double get_raw_left_velocity();
  double get_raw_right_velocity();
  std::vector<double> get_distance_readings();

  double get_actual_voltage();
  double get_current_draw();

  /**
   * @brief High-fidelity EKF estimates (SI units)
   */
  double get_linear_velocity();  // m/s
  double get_angular_velocity(); // rad/s (CW positive)

  double get_wheel_diameter() const { return wheelDiameter_m * 39.37; }
  double get_track_width() const { return trackWidth_m * 39.37; }

  void update();
  void reset_sensors();

  /**
   * @brief Physics-based dynamics prediction using blasfeo.
   */
  void predict_dynamics(double dt);

  /**
   * @brief High-performance STSMC velocity step.
   * Directly sets motor voltages.
   */
  void velocity_control(double v_ref, double w_ref, double dt);
  void set_gains(const GainsFG &gains);
  void set_stsmc_gains(const STSMCParams &params) { st_g_ = params; }
  void set_ekf_config(const EKFConfig &config) { ekf_config_ = config; }

private:
  std::vector<std::unique_ptr<pros::v5::Motor>> left_motors;
  std::vector<std::unique_ptr<pros::v5::Motor>> right_motors;
  std::vector<std::unique_ptr<pros::v5::Imu>> imus;
  std::vector<std::unique_ptr<pros::v5::Distance>> distances;

  std::vector<MotorData> left_data;
  std::vector<MotorData> right_data;
  std::vector<double> last_distances;

  EKFConfig ekf_config_;
  double P_cov[4] = {1.0, 0.0, 0.0, 1.0}; // 2x2 covariance matrix
  double gear_ratio;
  uint32_t last_time = 0;

  double last_left_pos = 0;
  double last_right_pos = 0;
  double last_left_vel = 0;
  double last_right_vel = 0;
  double last_heading = 0;

  double wheelDiameter_m;
  double trackWidth_m;

  GainsFG g_ = {0, 0, 0, 0, 0, 0, 0, 0};
  STSMCParams st_g_;

  // BLASFEO structures
  struct blasfeo_dmat sF, sG;
  struct blasfeo_dvec sx, su, sxdot;
  void *blasfeo_mem = nullptr;

  // STSMC State
  double z_int[2] = {0, 0};
  double x_ref_prev[2] = {0, 0};
  bool first_ctrl_step = true;

  void init_blasfeo();
};
