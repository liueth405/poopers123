#include "motions/motion.hpp"
#include <cmath>

MotionManager::MotionManager(Chassis &chassis, ParticleFilter &pf,
                             MPCC_Controller &mpcc_controller)
    : chassis(chassis), pf(pf), mpcc_controller(mpcc_controller) {}

void MotionManager::followPath(PathTarget target) {
  cancelMotion();
  current_path = target; // Store original (inches) for distance checking

  // Convert Waypoints to Meters for the MPCC Solver
  std::vector<PathVec2> waypoints_m;
  for (const auto &wp : target.waypoints) {
    waypoints_m.push_back({wp.x * 0.0254, wp.y * 0.0254});
  }
  mpcc_controller.setPath(waypoints_m);

  mpcc_controller.setControlBounds(target.mpc_control_bounds);

  // target.v_target_end is given in in/s, convert to m/s
  mpcc_controller.setTerminalTargetVelocity(target.v_target_end * 0.0254);
  mpcc_controller.setIntakeOffset(target.intake_offset *
                                  0.0254); // Convert inches to meters
  mpcc_controller.setWeights(target.weights);

  Estimate est = pf.estimate();
  MPCC_Controller::State start_state;
  start_state.x = est.x * 0.0254;
  start_state.y = est.y * 0.0254;
  start_state.theta = est.theta_deg * DEG_TO_RAD;
  start_state.vy = est.lateral * 0.0254;

  // Force start at the beginning of the path (s=0)
  start_state.s = 0.0;

  current_mpcc_s = start_state.s;
  prev_mpcc_u = {chassis.get_linear_velocity() * 0.0254,
                 chassis.get_angular_velocity(), 0.0};

  mpcc_controller.resetBuffer(start_state);

  state = MotionState::FOLLOWING_PATH;

  if (target.log_mpc_to_sd) {
    mpc_sd_log_.begin(target.mpc_log_csv_path);
  }
}

void MotionManager::waitUntilPathDone() {
  while (state == MotionState::FOLLOWING_PATH) {
    pros::delay(10);
  }
}

void MotionManager::moveTo(PointTarget target) {
  cancelMotion();
  current_point = target;
  linear_controller = PID(target.linear_pid);
  linear_controller.set_integral_limit(target.integral_limit);
  v_slew.limit = target.slew;
  v_slew.reset(0);
  state = MotionState::DRIVING;
}

void MotionManager::turnTo(TurnTarget target) {
  cancelMotion();
  current_turn = target;
  angular_controller = PID(target.angular_pid);
  angular_controller.set_integral_limit(target.integral_limit);
  w_slew.limit = target.slew;
  w_slew.reset(0);
  state = MotionState::TURNING;
}

void MotionManager::turnToPoint(double x, double y, TurnTarget::Direction dir,
                                PIDCoeffs angular_pid, double max_v,
                                double slew, double integral_limit) {
  Estimate est = pf.estimate();
  double dx = x - est.x;
  double dy = y - est.y;
  double target_angle = std::atan2(dx, dy) * RAD_TO_DEG;

  TurnTarget target;
  target.target_theta = target_angle;
  target.direction = dir;
  target.angular_pid = angular_pid;
  target.max_v = max_v;
  target.slew = slew;
  target.integral_limit = integral_limit;
  turnTo(target);
}

void MotionManager::cancelMotion() {
  mpc_sd_log_.end();
  state = MotionState::IDLE;
  chassis.tank(0, 0);
}

void MotionManager::waitUntilPoint(double x, double y, double margin) {
  while (!isDone()) {
    Estimate est = pf.estimate();
    double dist = std::sqrt(std::pow(x - est.x, 2) + std::pow(y - est.y, 2));
    if (dist < margin)
      break;
    pros::delay(10);
  }
}

void MotionManager::waitUntilAngle(double theta, double margin) {
  while (!isDone()) {
    Estimate est = pf.estimate();
    double error = std::abs(wrap_angle(theta - est.theta_deg));
    if (error < margin)
      break;
    pros::delay(10);
  }
}

void MotionManager::waitUntilPose(double x, double y, double theta,
                                  double pos_margin, double theta_margin) {
  while (!isDone()) {
    Estimate est = pf.estimate();
    double dist = std::sqrt(std::pow(x - est.x, 2) + std::pow(y - est.y, 2));
    double angle_error = std::abs(wrap_angle(theta - est.theta_deg));
    if (dist < pos_margin && angle_error < theta_margin)
      break;
    pros::delay(10);
  }
}

void MotionManager::waitUntilDistance(double dist) {
  Estimate start_est = pf.estimate();
  while (!isDone()) {
    Estimate curr = pf.estimate();
    double d = std::sqrt(std::pow(curr.x - start_est.x, 2) +
                         std::pow(curr.y - start_est.y, 2));
    if (d >= dist)
      break;
    pros::delay(10);
  }
}

void MotionManager::waitUntilLinearVelocity(double velocity, bool greaterThan) {
  while (!isDone()) {
    double current_v = std::abs(
        (chassis.get_left_velocity() + chassis.get_right_velocity()) / 2.0);
    if (greaterThan ? (current_v >= velocity) : (current_v <= velocity))
      break;
    pros::delay(10);
  }
}

void MotionManager::waitUntilAngularVelocity(double omega, bool greaterThan) {
  while (!isDone()) {
    double current_w = std::abs(chassis.get_angular_velocity());
    if (greaterThan ? (current_w >= omega) : (current_w <= omega))
      break;
    pros::delay(10);
  }
}

double MotionManager::wrap_angle(double angle) const {
  angle = std::fmod(angle + 180, 360);
  if (angle < 0)
    angle += 360;
  return angle - 180;
}

void MotionManager::loop() {
  uint32_t last_loop_time = pros::millis();

  uint32_t last_pf_time = pros::millis();
  uint32_t last_mpc_time = pros::millis();
  uint32_t last_vel_time = pros::millis();

  while (true) {
    chassis.update();

    static float prev_heading = 0;
    float curr_heading = chassis.get_rotation() * DEG_TO_RAD;
    float vl = chassis.get_left_velocity();
    float vr = chassis.get_right_velocity();
    std::vector<double> dists = chassis.get_distance_readings();
    std::vector<float> f_dists(dists.begin(), dists.end());

    uint32_t now_pf = pros::millis();
    double dt_pf = (now_pf - last_pf_time) / 1000.0;
    last_pf_time = now_pf;
    if (dt_pf <= 0.0 || dt_pf > 0.1)
      dt_pf = 0.01; // Safety bounds

    pf.update(curr_heading, prev_heading, vl, vr, f_dists.data(), dt_pf);
    prev_heading = curr_heading;
    Estimate est = pf.estimate();

    // Prepare target outputs for the single cascade call at the end
    double target_v = 0.0;
    double target_w = 0.0;

    if (state == MotionState::IDLE) {
      target_v = 0.0;
      target_w = 0.0;

    } else if (state == MotionState::DRIVING) {
      double dx = current_point.x - est.x;
      double dy = current_point.y - est.y;
      double dist = std::abs(std::sqrt(dx * dx + dy * dy));

      if (dist < current_point.tolerance) {
        state = MotionState::IDLE;
      } else {
        float angle_to_target = std::atan2(dx, dy) * RAD_TO_DEG;
        float angle_error = wrap_angle(angle_to_target - est.theta_deg);

        if (current_point.reversed) {
          angle_error = wrap_angle(angle_error + 180);
        }

        double error_rad = angle_error * DEG_TO_RAD;
        double K = 0;
        if (dist > current_point.exit_turn_dist) {
          K = (2.0 * std::sin(error_rad)) / std::max(dist, 2.0);
        }

        double v_req = linear_controller.update(dist, 0);
        v_req = std::clamp(v_req, -current_point.max_v, current_point.max_v);
        if (std::abs(v_req) < current_point.min_v)
          v_req = (v_req > 0 ? 1 : -1) * current_point.min_v;

        double v_slewed;
        if (std::abs(v_req) > std::abs(v_slew.current))
          v_slewed = v_slew.update(v_req);
        else
          v_slewed = v_slew.current = v_req;

        if (current_point.reversed)
          v_slewed *= -1;

        double w = std::abs(v_slewed) * K * current_point.steering_gain;
        double v_ms = v_slewed * 0.0254; // Convert in/s to m/s

        target_v = v_ms;
        target_w = w;
      }

    } else if (state == MotionState::TURNING) {
      double angle_error =
          wrap_angle(current_turn.target_theta - est.theta_deg);

      if (std::abs(angle_error) < current_turn.tolerance) {
        state = MotionState::IDLE;
      } else {
        if (std::abs(angle_error) > 15.0) {
          if (current_turn.direction == TurnTarget::Direction::LEFT &&
              angle_error > 0) {
            angle_error -= 360;
          } else if (current_turn.direction == TurnTarget::Direction::RIGHT &&
                     angle_error < 0) {
            angle_error += 360;
          }
        }

        double w_req = angular_controller.update(angle_error, 0);
        w_req = std::clamp(w_req, -current_turn.max_v, current_turn.max_v);
        if (std::abs(w_req) < current_turn.min_v)
          w_req = (w_req > 0 ? 1 : -1) * current_turn.min_v;

        double w_slewed;
        if (std::abs(w_req) > std::abs(w_slew.current))
          w_slewed = w_slew.update(w_req);
        else
          w_slewed = w_slew.current = w_req;

        double w_rads = w_slewed * DEG_TO_RAD;

        target_v = 0.0;
        target_w = w_rads;
      }

    } else if (state == MotionState::FOLLOWING_PATH) {
      double path_length_m = mpcc_controller.getPathLength();
      bool at_path_end =
          (current_mpcc_s >= path_length_m - 0.02); // 2cm tolerance
      auto last_wp_inches = current_path.waypoints.back();
      double real_dx = last_wp_inches.x - est.x;
      double real_dy = last_wp_inches.y - est.y;
      double real_dist_to_end =
          std::sqrt(real_dx * real_dx + real_dy * real_dy);

      // Calculate final orientation error (Wrapped correctly)
      Pose final_ref = mpcc_controller.getPathPoseAtS(path_length_m);
      double cur_th_rad = (est.theta_deg * DEG_TO_RAD);
      double eth_rad = cur_th_rad - final_ref.theta;
      while (eth_rad > M_PI)
        eth_rad -= 2.0 * M_PI;
      while (eth_rad < -M_PI)
        eth_rad += 2.0 * M_PI;
      double eth_deg = std::abs(eth_rad * (180.0 / M_PI));

      if ((real_dist_to_end < 2.0 && eth_deg < 5.0) ||
          (current_mpcc_s >= path_length_m - 0.01)) {
        mpc_sd_log_.end();
        state = MotionState::IDLE;
        chassis.tank(0, 0);
      } else {
        uint32_t now_mpc = pros::millis();
        double dt_mpc = (now_mpc - last_mpc_time) / 1000.0;
        last_mpc_time = now_mpc;
        if (dt_mpc < 1e-4)
          dt_mpc = 0.01;
        if (dt_mpc > 0.1)
          dt_mpc = 0.01;

        Estimate est = pf.estimate();
        double v_cur = chassis.get_linear_velocity() * 0.0254; // m/s
        double w_cur = chassis.get_angular_velocity();
        double dt_latency = dt_mpc;

        MPCC_Controller::State mpc_state;
        mpc_state.theta = (est.theta_deg * DEG_TO_RAD) + (w_cur * dt_latency);
        mpc_state.x =
            (est.x * 0.0254) + (v_cur * std::sin(mpc_state.theta) * dt_latency);
        mpc_state.y =
            (est.y * 0.0254) + (v_cur * std::cos(mpc_state.theta) * dt_latency);
        mpc_state.vy = est.lateral * 0.0254;
        mpc_state.s = current_mpcc_s;

        double s_predicted = current_mpcc_s + prev_mpcc_u.vs * dt_mpc;
        s_predicted =
            std::clamp(s_predicted, 0.0, mpcc_controller.getPathLength());
        double offset_m = current_path.intake_offset * 0.0254;
        double xi = mpc_state.x + offset_m * std::sin(mpc_state.theta);
        double yi = mpc_state.y + offset_m * std::cos(mpc_state.theta);
        double s_proj_raw =
            mpcc_controller.projectSFromPosition(xi, yi, current_mpcc_s);
        double s_proj = std::max(current_mpcc_s, s_proj_raw);

        double dist_to_path_end = path_length_m - current_mpcc_s;
        double progress_alpha = (dist_to_path_end < 0.2) ? 0.05 : 0.2;
        double s_fused =
            progress_alpha * s_proj + (1.0 - progress_alpha) * s_predicted;

        double max_lead = 0.20; // 8 inches max pull ahead
        s_fused = std::min(s_predicted, s_proj + max_lead);

        current_mpcc_s = std::max(current_mpcc_s, s_fused);
        mpc_state.s = current_mpcc_s;

        mpcc_controller.updateParameters(prev_mpcc_u, mpc_state.s, dt_mpc);

        // double end_scale = std::clamp(dist_to_path_end / 0.3, 0.4,
        //                               1.0); // Don't scale below 40%
        // MPCC_Controller::ControlBounds adjusted_bounds =
        //     current_path.mpc_control_bounds;

        // adjusted_bounds.v_max *= end_scale;
        // adjusted_bounds.vs_max *= end_scale;

        // mpcc_controller.setControlBounds(adjusted_bounds);
        mpcc_controller.setControlBounds(current_path.mpc_control_bounds);

        if (dist_to_path_end < 0.04) { // (~1.5 inches)
          mpc_sd_log_.end();
          state = MotionState::IDLE;
          chassis.tank(0, 0);
          return;
        }

        prev_mpcc_u = mpcc_controller.runStep(mpc_state, prev_mpcc_u);
        MPCC_Controller::Control u = prev_mpcc_u;

        if (!std::isfinite(u.v) || !std::isfinite(u.w) ||
            !std::isfinite(u.vs)) {
          u = {0.0, 0.0, 0.0};
        }
        prev_mpcc_u = u;

        static MPCC_Controller::Control last_output = {0.0, 0.0, 0.0};
        u.w = std::clamp(u.w, last_output.w - current_path.max_dw_per_tick,
                         last_output.w + current_path.max_dw_per_tick);
        u.v = std::clamp(u.v, last_output.v - current_path.max_dv_per_tick,
                         last_output.v + current_path.max_dv_per_tick);
        last_output = u;
        target_v = u.v;
        target_w = u.w;

        if (mpc_sd_log_.is_open()) {
          const double meas_v_ins = chassis.get_linear_velocity();
          const double meas_w_rads = chassis.get_angular_velocity();
          const MPCC_Controller::SolveInfo si =
              mpcc_controller.getLastSolveInfo();
          const Pose refp = mpcc_controller.getPathPoseAtS(mpc_state.s);
          const Pose refp_proj = mpcc_controller.getPathPoseAtS(s_proj);
          mpc_sd_log_.log_row(
              pros::millis(), mpcc_controller.getPathLength(), est.x, est.y,
              est.theta_deg, est.lateral, mpc_state.x, mpc_state.y,
              mpc_state.theta, mpc_state.vy, mpc_state.s, s_proj, 0.0,
              current_mpcc_s, u.v, u.w, u.vs, si.status, si.time_ms,
              si.sqp_iter, si.kkt_norm, refp.x, refp.y, refp.theta, meas_v_ins,
              meas_w_rads, target_v, target_w);
        }
      }
    }

    uint32_t now_vel = pros::millis();
    double dt_vel = (now_vel - last_vel_time) / 1000.0;
    last_vel_time = now_vel;
    if (dt_vel <= 0.0 || dt_vel > 0.1)
      dt_vel = 0.01; // Safety bounds
    chassis.velocity_control(target_v, target_w, dt_vel);

    pros::Task::delay_until(&last_loop_time, 10);
  }
}