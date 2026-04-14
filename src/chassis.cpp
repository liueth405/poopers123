#include "chassis.hpp"
#include <cmath>
#include <cstdint>

Chassis::Chassis(std::vector<int> left_ports, std::vector<int> right_ports,
                 std::vector<int> imu_ports, std::vector<int> distance_ports,
                 double gear_ratio, double wheelDiameter, double trackWidth,
                 pros::v5::MotorGears cartridge, EKFConfig config)
    : gear_ratio(gear_ratio), wheelDiameter_m(wheelDiameter * 0.0254),
      trackWidth_m(trackWidth * 0.0254), ekf_config_(config) {

  first_ctrl_step = true;
  st_g_ = STSMCParams(); // Load defaults

  for (int port : left_ports) {
    auto m = std::make_unique<pros::v5::Motor>(
        port, cartridge, pros::v5::MotorEncoderUnits::degrees);
    left_motors.push_back(std::move(m));
    left_data.push_back(MotorData());
  }

  for (int port : right_ports) {
    auto m = std::make_unique<pros::v5::Motor>(
        port, cartridge, pros::v5::MotorEncoderUnits::degrees);
    right_motors.push_back(std::move(m));
    right_data.push_back(MotorData());
  }

  for (int port : imu_ports) {
    imus.push_back(std::make_unique<pros::v5::Imu>(port));
  }

  for (int port : distance_ports) {
    distances.push_back(
        std::make_unique<V5_DeviceT>(vexDeviceGetByIndex(port - 1)));
    last_times_distance.push_back(0);
    last_distances.push_back(0);
  }

  init_blasfeo();
  last_time = pros::millis();
}

Chassis::~Chassis() {
  if (blasfeo_mem)
    free(blasfeo_mem);
}

void Chassis::init_blasfeo() {
  int m = 2, n = 2;
  size_t size = 0;
  size += blasfeo_memsize_dmat(m, n);  // sF
  size += blasfeo_memsize_dmat(m, n);  // sG
  size += blasfeo_memsize_dvec(m) * 3; // sx, su, sxdot
  size += 128;                         // alignment padding

  blasfeo_mem = malloc(size);
  char *ptr = (char *)blasfeo_mem;

  blasfeo_create_dmat(m, n, &sF, ptr);
  ptr += blasfeo_memsize_dmat(m, n);

  blasfeo_create_dmat(m, n, &sG, ptr);
  ptr += blasfeo_memsize_dmat(m, n);

  blasfeo_create_dvec(m, &sx, ptr);
  ptr += blasfeo_memsize_dvec(m);

  blasfeo_create_dvec(m, &su, ptr);
  ptr += blasfeo_memsize_dvec(m);

  blasfeo_create_dvec(m, &sxdot, ptr);

  // Initialize matrices to zero
  blasfeo_dgese(m, n, 0.0, &sF, 0, 0);
  blasfeo_dgese(m, n, 0.0, &sG, 0, 0);
  blasfeo_dvecse(m, 0.0, &sx, 0);
  blasfeo_dvecse(m, 0.0, &su, 0);
}

void Chassis::set_gains(const GainsFG &gains) {
  g_ = gains;
  // Pack into blasfeo mats
  // F matrix (dx/dt = Fx + Gu)
  // F = [-lambda_v, 0; 0, -lambda_w]
  // Packing Column-Major: [R0C0, R1C0, R0C1, R1C1]
  double f_data[4] = {-(g_.lambda_v), 0, 0, -(g_.lambda_w)};
  blasfeo_pack_dmat(2, 2, f_data, 2, &sF, 0, 0);

  // G matrix
  // x_dot = Fx + Gu
  // u = [uL, uR]
  // v_dot = alpha_v * 0.5 * (uL + uR) - lambda_v * v
  // w_dot = alpha_w * 0.5 * (uL - uR) - lambda_w * w
  // Packing Column-Major: [R0C0, R1C0, R0C1, R1C1]
  double g_data[4] = {g_.alpha_v * 0.5, g_.alpha_w * 0.5, g_.alpha_v * 0.5,
                      -g_.alpha_w * 0.5};
  blasfeo_pack_dmat(2, 2, g_data, 2, &sG, 0, 0);
}

void Chassis::tank(int left_voltage, int right_voltage) {
  for (size_t i = 0; i < left_motors.size(); ++i) {
    left_motors[i]->move_voltage(left_voltage);
    left_data[i].voltage = (double)left_voltage;
  }
  for (size_t i = 0; i < right_motors.size(); ++i) {
    right_motors[i]->move_voltage(right_voltage);
    right_data[i].voltage = (double)right_voltage;
  }
}

void Chassis::update() {
  uint32_t now = pros::millis();
  double dt = (now - last_time) / 1000.0;
  if (dt <= 0)
    return;
  last_time = now;

  // 1. Get Measurements and Commands
  double vlt_l = 0, vlt_r = 0;
  int cl_v = 0, cr_v = 0;
  for (auto &m : left_motors) {
    double v = m->get_voltage();
    if (v != PROS_ERR_F) {
      vlt_l += v;
      cl_v++;
    }
  }
  if (cl_v > 0)
    vlt_l /= (cl_v * 1000.0);

  for (auto &m : right_motors) {
    double v = m->get_voltage();
    if (v != PROS_ERR_F) {
      vlt_r += v;
      cr_v++;
    }
  }
  if (cr_v > 0)
    vlt_r /= (cr_v * 1000.0);

  double vl_meas = 0, vr_meas = 0;
  int cl = 0, cr = 0;

  for (auto &m : left_motors) {
    double v = m->get_actual_velocity();
    if (v != PROS_ERR_F) {
      vl_meas += v;
      cl++;
    }
  }
  if (cl > 0)
    vl_meas /= cl;

  for (auto &m : right_motors) {
    double v = m->get_actual_velocity();
    if (v != PROS_ERR_F) {
      vr_meas += v;
      cr++;
    }
  }
  if (cr > 0)
    vr_meas /= cr;

  // Convert RPM -> m/s
  double vl_si = (vl_meas / gear_ratio) * (M_PI / 60.0) * wheelDiameter_m;
  double vr_si = (vr_meas / gear_ratio) * (M_PI / 60.0) * wheelDiameter_m;

  double v_enc = (vl_si + vr_si) / 2.0;          // m/s
  double w_enc = (vl_si - vr_si) / trackWidth_m; // rad/s (CW positive)

  double prev_heading = last_heading;
  double curr_heading = get_rotation();
  last_heading = curr_heading;

  double dh = curr_heading - prev_heading;
  while (dh > 180)
    dh -= 360;
  while (dh < -180)
    dh += 360;

  double w_imu = dh * (M_PI / 180.0) / dt;

  // 2. Physics Prediction Step
  double u_data[2] = {vlt_l, vlt_r};
  blasfeo_pack_dvec(2, u_data, 1, &su, 0);

  // Get previous states
  double v_prev = blasfeo_dvecex1(&sx, 0);
  double w_prev = blasfeo_dvecex1(&sx, 1);

  // Calculate dynamic prediction (dx/dt = F*sx + G*su)
  blasfeo_dgemv_n(2, 2, 1.0, &sF, 0, 0, &sx, 0, 0.0, &sxdot, 0, &sxdot, 0);
  blasfeo_dgemv_n(2, 2, 1.0, &sG, 0, 0, &su, 0, 1.0, &sxdot, 0, &sxdot, 0);

  // =========================================================================
  // NEW: ROBUST SLIP & STALL PROTECTION
  // =========================================================================
  const double MAX_DV =
      (g_.max_accel_wall_v + g_.max_accel_wall_w * trackWidth_m / 2.0) * dt;
  bool is_slipping_or_stalled = false;

  // A. Wheel Spin-out Protection: Clamp impossible encoder accelerations
  if (v_enc > v_prev + MAX_DV) {
    v_enc = v_prev + MAX_DV;
    is_slipping_or_stalled = true;
  } else if (v_enc < v_prev - MAX_DV) {
    v_enc = v_prev - MAX_DV;
    is_slipping_or_stalled = true;
  }

  // B. Stall / Climbing Protection: Encoder is slow, but voltage predicts fast
  double expected_v_from_voltage = v_prev + blasfeo_dvecex1(&sxdot, 0) * dt;
  if (std::abs(v_enc) < 0.15 && std::abs(expected_v_from_voltage) > 0.5) {
    is_slipping_or_stalled = true;
  }

  double v_pred, w_pred;

  if (is_slipping_or_stalled) {
    // "Just use the previous velocity": Bypass the voltage model completely
    // so the prediction doesn't push the robot through a wall mathematically.
    v_pred = v_prev * (1.0 - g_.lambda_v * dt);
    w_pred = w_prev * (1.0 - g_.lambda_w * dt);

    // Sync the blasfeo vector to match this override
    blasfeo_dvecin1(v_pred, &sx, 0);
    blasfeo_dvecin1(w_pred, &sx, 1);
  } else {
    // Normal Prediction Step
    blasfeo_dvecad(2, dt, &sxdot, 0, &sx, 0);
    v_pred = blasfeo_dvecex1(&sx, 0);
    w_pred = blasfeo_dvecex1(&sx, 1);
  }
  // =========================================================================

  // EKF Covariance Prediction
  double a00 = 1.0 - g_.lambda_v * dt;
  double a11 = 1.0 - g_.lambda_w * dt;
  P_cov[0] = a00 * a00 * P_cov[0] + ekf_config_.Q_v * dt;
  P_cov[1] = a00 * a11 * P_cov[1];
  P_cov[2] = a00 * a11 * P_cov[2];
  P_cov[3] = a11 * a11 * P_cov[3] + ekf_config_.Q_w * dt;

  // EKF Measurement Updates (REMOVED MAHALANOBIS REJECTION)
  // Since we already clamped the encoders to physical reality, we ALWAYS
  // want the EKF to trust them now.
  auto apply_update = [&](double z, double h0, double h1, double R) {
    double z_hat = h0 * v_pred + h1 * w_pred;
    double y = z - z_hat;

    double S = h0 * h0 * P_cov[0] + h0 * h1 * P_cov[1] + h1 * h0 * P_cov[2] +
               h1 * h1 * P_cov[3] + R;

    double K0 = (P_cov[0] * h0 + P_cov[1] * h1) / S;
    double K1 = (P_cov[2] * h0 + P_cov[3] * h1) / S;

    v_pred += K0 * y;
    w_pred += K1 * y;

    double t00 = 1.0 - K0 * h0, t01 = -K0 * h1;
    double t10 = -K1 * h0, t11 = 1.0 - K1 * h1;

    double p0_new = t00 * P_cov[0] + t01 * P_cov[1];
    double p1_new = t10 * P_cov[0] + t11 * P_cov[1];
    double p2_new = t00 * P_cov[2] + t01 * P_cov[3];
    double p3_new = t10 * P_cov[2] + t11 * P_cov[3];

    P_cov[0] = p0_new;
    P_cov[1] = p1_new;
    P_cov[2] = p2_new;
    P_cov[3] = p3_new;
  };

  // Process measurements sequentially
  apply_update(w_imu, 0.0, 1.0, ekf_config_.R_w_imu);
  apply_update(w_enc, 0.0, 1.0, ekf_config_.R_w_enc);
  apply_update(v_enc, 1.0, 0.0, ekf_config_.R_v_enc);

  // Write updated states back to blasfeo vector sx
  blasfeo_dvecin1(v_pred, &sx, 0);
  blasfeo_dvecin1(w_pred, &sx, 1);

  double v_new = v_pred;
  double w_new = w_pred;

  // 4. Update External Outputs
  double vl_final = v_new + w_new * (trackWidth_m / 2.0);
  double vr_final = v_new - w_new * (trackWidth_m / 2.0);

  last_left_vel = vl_final / 0.0254;  // -> in/s
  last_right_vel = vr_final / 0.0254; // -> in/s

  // Position: use raw position (scaled)
  double pl_avg = 0;
  int plc = 0;
  for (auto &m : left_motors) {
    double p = m->get_position();
    if (p != PROS_ERR_F) {
      pl_avg += p;
      plc++;
    }
  }
  if (plc > 0)
    last_left_pos = ((pl_avg / plc) / gear_ratio) * (M_PI * wheelDiameter_m) /
                    360.0 / 0.0254;

  double pr_avg = 0;
  int prc = 0;
  for (auto &m : right_motors) {
    double p = m->get_position();
    if (p != PROS_ERR_F) {
      pr_avg += p;
      prc++;
    }
  }
  if (prc > 0)
    last_right_pos = ((pr_avg / prc) / gear_ratio) * (M_PI * wheelDiameter_m) /
                     360.0 / 0.0254;

  for (size_t i = 0; i < distances.size(); ++i) {
    uint32_t time = vexDeviceGetTimestamp(*distances[i]);
    if (time == last_times_distance[i]) {
      last_distances[i] = 9999;
      continue;
    }
    last_times_distance[i] = time;
    double d = (double)vexDeviceDistanceDistanceGet(*distances[i]) / 25.4;
    if (d != PROS_ERR_F)
      last_distances[i] = d;
  }
}

double Chassis::get_rotation() {
  double heading_sum = 0;
  int connected_count = 0;
  for (auto &imu : imus) {
    double cur = imu->get_rotation();
    if (imu->is_installed() && cur != PROS_ERR_F && !std::isnan(cur) &&
        imu->get_status() != pros::ImuStatus::calibrating) {
      heading_sum += cur;
      connected_count++;
    }
  }
  if (connected_count > 0)
    last_heading = heading_sum / connected_count; // CW positive
  return last_heading;
}

double Chassis::get_left_position() { return last_left_pos; }
double Chassis::get_right_position() { return last_right_pos; }
double Chassis::get_left_velocity() {
  double v = blasfeo_dvecex1(&sx, 0);
  double w = blasfeo_dvecex1(&sx, 1);
  return (v + w * (trackWidth_m / 2.0)) * 39.3701; // Convert m/s -> in/s
}

double Chassis::get_right_velocity() {
  double v = blasfeo_dvecex1(&sx, 0);
  double w = blasfeo_dvecex1(&sx, 1);
  return (v - w * (trackWidth_m / 2.0)) * 39.3701; // Convert m/s -> in/s
}

double Chassis::get_raw_left_velocity() {
  double vl_meas = 0;
  int cl = 0;
  for (auto &m : left_motors) {
    double v = m->get_actual_velocity();
    if (v != PROS_ERR_F) {
      vl_meas += v;
      cl++;
    }
  }
  if (cl > 0)
    vl_meas /= cl;

  double vl_si = (vl_meas / gear_ratio) * (M_PI / 60.0) * wheelDiameter_m;
  return vl_si / 0.0254; // in/s
}

double Chassis::get_raw_right_velocity() {
  double vr_meas = 0;
  int cr = 0;
  for (auto &m : right_motors) {
    double v = m->get_actual_velocity();
    if (v != PROS_ERR_F) {
      vr_meas += v;
      cr++;
    }
  }
  if (cr > 0)
    vr_meas /= cr;

  double vr_si = (vr_meas / gear_ratio) * (M_PI / 60.0) * wheelDiameter_m;
  return vr_si / 0.0254; // in/s
}

double Chassis::get_actual_voltage() {
  double sum = 0;
  int count = 0;
  for (auto &m : left_motors) {
    double v = m->get_voltage();
    if (v != PROS_ERR_F) {
      sum += std::abs(v);
      count++;
    }
  }
  for (auto &m : right_motors) {
    double v = m->get_voltage();
    if (v != PROS_ERR_F) {
      sum += std::abs(v);
      count++;
    }
  }
  return count > 0 ? (sum / count) / 1000.0 : 0;
}

double Chassis::get_current_draw() {
  double sum = 0;
  int count = 0;
  for (auto &m : left_motors) {
    double c = m->get_current_draw();
    if (c != PROS_ERR_F) {
      sum += std::abs(c);
      count++;
    }
  }
  for (auto &m : right_motors) {
    double c = m->get_current_draw();
    if (c != PROS_ERR_F) {
      sum += std::abs(c);
      count++;
    }
  }
  return count > 0 ? (sum / count) / 1000.0 : 0;
}

double Chassis::get_linear_velocity() {
  return blasfeo_dvecex1(&sx, 0) * 39.3701; // in/s
}

double Chassis::get_angular_velocity() {
  return blasfeo_dvecex1(&sx, 1); // rad/s
}

std::vector<double> Chassis::get_distance_readings() {
  std::vector<double> readings;
  for (size_t i = 0; i < distances.size(); ++i)
    readings.push_back(last_distances[i]); // already in inches
  return readings;
}
void Chassis::reset_sensors() {
  for (auto &m : left_motors)
    m->tare_position();
  for (auto &m : right_motors)
    m->tare_position();
  for (auto &imu : imus)
    imu->reset(true);
  blasfeo_dvecse(2, 0.0, &sx, 0);
  last_left_pos = 0;
  last_right_pos = 0;
  last_left_vel = 0;
  last_right_vel = 0;
  last_heading = 0;
}

void Chassis::velocity_control(double v_ref, double w_ref, double dt) {

  // 0. Handle initial state
  if (first_ctrl_step) {
    x_ref_prev[0] = v_ref;
    x_ref_prev[1] = w_ref;
    first_ctrl_step = false;
  }

  // 1. Reference Derivative (Feedforward)
  // The slew rate in MotionManager ensures v_ref is smooth, so differentiation
  // is safe!
  double v_ref_dot = (v_ref - x_ref_prev[0]) / dt;
  double w_ref_dot = (w_ref - x_ref_prev[1]) / dt;

  x_ref_prev[0] = v_ref;
  x_ref_prev[1] = w_ref;

  // 2. Tracking Error and Sliding Surface
  // v_meas = blasfeo_dvecex1(&sx, 0), w_meas = blasfeo_dvecex1(&sx, 1)
  double e_v = blasfeo_dvecex1(&sx, 0) - v_ref;
  double e_w = blasfeo_dvecex1(&sx, 1) - w_ref;

  double s_v = st_g_.Lambda[0] * e_v;
  double s_w = st_g_.Lambda[1] * e_w;

  // 3. Super-Twisting Control
  auto smooth_sgn = [&](double s) {
    if (s > st_g_.s_eps)
      return 1.0;
    if (s < -st_g_.s_eps)
      return -1.0;
    return s / st_g_.s_eps;
  };

  double sgn_v = smooth_sgn(s_v);
  double sgn_w = smooth_sgn(s_w);

  // Integrators
  // Integrators with "Leakage" for 0-target stopping
  // If target is 0, decay the integrator by a factor of 5 (bleeds to 0 in ~0.5
  // seconds)
  double leak_v = (std::abs(v_ref) < 0.01) ? 5.0 * z_int[0] : 0.0;
  z_int[0] += (-st_g_.K2[0] * sgn_v - leak_v) * dt;
  z_int[0] = std::clamp(z_int[0], -st_g_.z_max[0], st_g_.z_max[0]);

  double leak_w = (std::abs(w_ref) < 0.01) ? 5.0 * z_int[1] : 0.0;
  z_int[1] += (-st_g_.K2[1] * sgn_w - leak_w) * dt;
  z_int[1] = std::clamp(z_int[1], -st_g_.z_max[1], st_g_.z_max[1]);

  // u_st = -K1 * sqrt(|s|) * sign(s) + z
  double u_st_v = -st_g_.K1[0] * std::sqrt(std::abs(s_v)) * sgn_v + z_int[0];
  double u_st_w = -st_g_.K1[1] * std::sqrt(std::abs(s_w)) * sgn_w + z_int[1];

  // 4. Model Decoupling (Feedback Linearization)
  // FIX: Use v_ref and w_ref instead of measurements to provide instant
  // steady-state Feedforward. This completely eliminates SMC Integrator windup
  // on step inputs.
  double xdot_des_v = v_ref_dot + g_.lambda_v * v_ref + u_st_v;
  double xdot_des_w = w_ref_dot + g_.lambda_w * w_ref + u_st_w;

  // 5. Inverse Dynamics Mapping
  double uL = (xdot_des_v / g_.alpha_v) + (xdot_des_w / g_.alpha_w);
  double uR = (xdot_des_v / g_.alpha_v) - (xdot_des_w / g_.alpha_w);

  // 6. Friction Feedforward
  // Calculate the physical target velocity of each wheel
  double vL_ref = v_ref + w_ref * (trackWidth_m / 2.0);
  double vR_ref = v_ref - w_ref * (trackWidth_m / 2.0);

  // a) Apply basic Coulomb friction to each wheel based on its rolling
  // direction
  if (std::abs(vL_ref) > 0.001)
    uL += (vL_ref > 0 ? g_.kS_v : -g_.kS_v);
  if (std::abs(vR_ref) > 0.001)
    uR += (vR_ref > 0 ? g_.kS_v : -g_.kS_v);

  // b) Apply Curvature-Scaled Scrub Friction
  // Calculate how much of the robot's movement is turning vs driving straight
  double turn_comp = std::abs(w_ref * (trackWidth_m / 2.0));
  double drive_comp = std::abs(v_ref);
  double turn_ratio = turn_comp / (drive_comp + turn_comp + 1e-6);

  // Scale the scrub friction by the turn ratio (wide arcs have less scrub than
  // point turns)
  double scrub_ff = std::max(0.0, g_.kS_w - g_.kS_v);
  if (std::abs(w_ref) > 0.001) {
    double scrub_volts = scrub_ff * turn_ratio;
    uL += (w_ref > 0 ? scrub_volts : -scrub_volts);
    uR += (w_ref > 0 ? -scrub_volts : scrub_volts);
  }

  // 7. Saturation (Uniform Scaling to preserve turn curvature)
  double max_vlt = std::max(std::abs(uL), std::abs(uR));
  if (max_vlt > 12.0) {
    double scale = 12.0 / max_vlt;
    uL *= scale;
    uR *= scale;
  }

  tank((int)(uL * 1000.0), (int)(uR * 1000.0));
}
