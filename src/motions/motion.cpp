#include "motions/motion.hpp"

MotionManager::MotionManager(Chassis &chassis, ParticleFilter &pf)
    : chassis(chassis), pf(pf) {}

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

  while (true) {
    chassis.update();

    static float prev_heading = 0;
    float curr_heading = chassis.get_rotation() * DEG_TO_RAD;
    float vl = chassis.get_left_velocity();
    float vr = chassis.get_right_velocity();
    std::vector<double> dists = chassis.get_distance_readings();
    std::vector<float> f_dists(dists.begin(), dists.end());

    pf.update(curr_heading, prev_heading, vl, vr, f_dists.data(), 0.01);
    prev_heading = curr_heading;
    Estimate est = pf.estimate();

    if (state == MotionState::IDLE) {
      chassis.velocity_control(0.0, 0.0, 0.01);
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
          // Clamp distance to 2 inches for the divisor to prevent steering explosion near target
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

        // Use absolute velocity for steering magnitude so the sign 
        // is determined solely by the (already adjusted) angle error.
        double w = std::abs(v_slewed) * K * current_point.steering_gain; // natural unit: rad/s

        // Cascade PID: Pass target velocities to the inner STSMC loop
        double v_ms = v_slewed * 0.0254; // Convert in/s to m/s
        chassis.velocity_control(v_ms, w, 0.01);
      }
    } else if (state == MotionState::TURNING) {
      double angle_error =
          wrap_angle(current_turn.target_theta - est.theta_deg);

      if (std::abs(angle_error) < current_turn.tolerance) {
        state = MotionState::IDLE;
      } else {
        // Only force direction if we are far away. 
        // When close (< 15 deg), always use shortest path to correct settlement.
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
          
        // Cascade PID: Convert deg/s to rad/s and pass to STSMC
        double w_rads = w_slewed * DEG_TO_RAD;
        chassis.velocity_control(0.0, w_rads, 0.01);
      }
    }

    pros::Task::delay_until(&last_loop_time, 10);
  }
}
