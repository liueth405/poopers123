#include "chassis.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>

Chassis::Chassis(std::vector<int> left_ports, std::vector<int> right_ports,
                 std::vector<int> imu_ports, std::vector<int> distance_ports,
                 double gear_ratio, pros::v5::MotorGears cartridge)
    : gear_ratio(gear_ratio) {

  for (int port : left_ports) {
    left_motors.push_back(std::make_unique<pros::v5::Motor>(port, cartridge));
    left_data.push_back(MotorData());
  }

  for (int port : right_ports) {
    right_motors.push_back(std::make_unique<pros::v5::Motor>(port, cartridge));
    right_data.push_back(MotorData());
  }

  for (int port : imu_ports) {
    imus.push_back(std::make_unique<pros::v5::Imu>(port));
  }

  for (int port : distance_ports) {
    distances.push_back(std::make_unique<pros::v5::Distance>(port));
    last_distances.push_back(0);
  }

  last_time = pros::millis();
}

void Chassis::tank(int left_voltage, int right_voltage) {
  for (auto &motor : left_motors) {
    motor->move_voltage(left_voltage);
  }
  for (auto &motor : right_motors) {
    motor->move_voltage(right_voltage);
  }
}

double Chassis::get_rotation() {
  double heading_sum = 0;
  int connected_count = 0;

  for (auto &imu : imus) {
    double current_heading = imu->get_rotation();
    if (imu->is_installed() && current_heading != PROS_ERR_F &&
        !std::isnan(current_heading)) {
      if (imu->get_status() != pros::ImuStatus::calibrating) {
        heading_sum += current_heading;
        connected_count++;
      }
    }
  }

  if (connected_count > 0) {
    last_heading = heading_sum / connected_count;
  }

  return last_heading;
}

double Chassis::get_left_position() {
  double sum = 0;
  int count = 0;
  for (const auto &data : left_data) {
    if (data.connected) {
      sum += data.position;
      count++;
    }
  }
  if (count > 0)
    last_left_pos = sum / count;
  return last_left_pos;
}

double Chassis::get_right_position() {
  double sum = 0;
  int count = 0;
  for (const auto &data : right_data) {
    if (data.connected) {
      sum += data.position;
      count++;
    }
  }
  if (count > 0)
    last_right_pos = sum / count;
  return last_right_pos;
}

double Chassis::get_left_velocity() {
  double sum = 0;
  int count = 0;
  for (const auto &data : left_data) {
    if (data.connected) {
      sum += data.velocity;
      count++;
    }
  }
  if (count > 0)
    last_left_vel = sum / count;
  return last_left_vel;
}

double Chassis::get_right_velocity() {
  double sum = 0;
  int count = 0;
  for (const auto &data : right_data) {
    if (data.connected) {
      sum += data.velocity;
      count++;
    }
  }
  if (count > 0)
    last_right_vel = sum / count;
  return last_right_vel;
}

std::vector<double> Chassis::get_distance_readings() {
  std::vector<double> readings;
  for (size_t i = 0; i < distances.size(); ++i) {
    double d = (double)distances[i]->get();
    if (d != PROS_ERR_F) {
      last_distances[i] = d;
    }
    readings.push_back(last_distances[i]);
  }
  return readings;
}

const Chassis::MotorData &Chassis::get_left_motor_data(size_t index) const {
  return left_data.at(index);
}

const Chassis::MotorData &Chassis::get_right_motor_data(size_t index) const {
  return right_data.at(index);
}

void Chassis::update() {
  uint32_t now = pros::millis();
  double dt = (now - last_time) / 1000.0;
  if (dt <= 0)
    return;
  last_time = now;

  get_rotation();

  auto update_motor_data = [&](auto &motors, auto &data_vec) {
    for (size_t i = 0; i < motors.size(); ++i) {
      double current_pos = motors[i]->get_position();

      if (current_pos != PROS_ERR_F) {
        if (!data_vec[i].connected) {
          data_vec[i].est_pos = current_pos;
          data_vec[i].est_vel = 0;
          data_vec[i].connected = true;
        }

        // alpha-beta filter implementation
        double pred_pos = data_vec[i].est_pos + data_vec[i].est_vel * dt;
        double pred_vel = data_vec[i].est_vel;

        // correction step
        double residual = current_pos - pred_pos;
        data_vec[i].est_pos = pred_pos + alpha * residual;
        data_vec[i].est_vel = pred_vel + (beta / dt) * residual;

        // convert deg/sec to RPM (RPM = (deg/sec) / 6.0)
        data_vec[i].velocity = (data_vec[i].est_vel / 6.0) * gear_ratio;
        data_vec[i].position = current_pos * gear_ratio;
      } else {
        data_vec[i].connected = false;
      }
    }
  };

  update_motor_data(left_motors, left_data);
  update_motor_data(right_motors, right_data);

  // Update distances
  for (size_t i = 0; i < distances.size(); ++i) {
    double d = (double)distances[i]->get();
    if (d != PROS_ERR_F) {
      last_distances[i] = d;
    }
  }
}

void Chassis::reset_sensors() {
  for (size_t i = 0; i < left_motors.size(); ++i) {
    left_motors[i]->tare_position();
    left_data[i].position = 0;
    left_data[i].velocity = 0;
    left_data[i].est_pos = 0;
    left_data[i].est_vel = 0;
  }
  for (size_t i = 0; i < right_motors.size(); ++i) {
    right_motors[i]->tare_position();
    right_data[i].position = 0;
    right_data[i].velocity = 0;
    right_data[i].est_pos = 0;
    right_data[i].est_vel = 0;
  }
  for (auto &imu : imus) {
    imu->tare_rotation();
  }
  for (size_t i = 0; i < last_distances.size(); ++i) {
    last_distances[i] = 0;
  }
  last_left_pos = 0;
  last_right_pos = 0;
  last_left_vel = 0;
  last_right_vel = 0;
  last_heading = 0;
}
