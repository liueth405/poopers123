#pragma once

#include "api.h"
#include <memory>
#include <vector>

/**
 * @brief Chassis class to handle robot movement and sensor data.
 */
class Chassis {
public:
  struct MotorData {
    double position = 0; // raw individual position (scaled by gear ratio)
    double velocity = 0; // filtered velocity estimate (RPM)

    // internal filter state (alpha-beta filter)
    double est_pos = 0;
    double est_vel = 0;
    bool connected = false;
  };

  /**
   * @brief Construct a new Chassis object.
   *
   * @param left_ports Ports of the left motors. Negative for reversed.
   * @param right_ports Ports of the right motors. Negative for reversed.
   * @param imu_ports Ports of the IMU(s).
   * @param distance_ports Ports of the Distance sensor(s).
   * @param gear_ratio External gear ratio (e.g., 0.6 for 36:60).
   * @param cartridge Motor cartridge (e.g., pros::v5::MotorGears::blue).
   */
  Chassis(std::vector<int> left_ports, std::vector<int> right_ports,
          std::vector<int> imu_ports, std::vector<int> distance_ports,
          double gear_ratio, pros::v5::MotorGears cartridge);

  /**
   * @brief Drive the robot using tank controls.
   *
   * @param left_voltage Voltage for left motors (-12000 to 12000).
   * @param right_voltage Voltage for right motors (-12000 to 12000).
   */
  void tank(int left_voltage, int right_voltage);

  /**
   * @brief Get the average rotation from all IMUs.
   */
  double get_rotation();

  /**
   * @brief Get the average position of the left motors.
   */
  double get_left_position();

  /**
   * @brief Get the average position of the right motors.
   */
  double get_right_position();

  /**
   * @brief Get the average filtered velocity of the left motors.
   */
  double get_left_velocity();

  /**
   * @brief Get the average filtered velocity of the right motors.
   */
  double get_right_velocity();

  /**
   * @brief Get distance sensor readings in mm.
   */
  std::vector<double> get_distance_readings();

  /**
   * @brief Get individual data for a left motor.
   */
  const MotorData &get_left_motor_data(size_t index) const;

  /**
   * @brief Get individual data for a right motor.
   */
  const MotorData &get_right_motor_data(size_t index) const;

  /**
   * @brief Updates the internal state, using an alpha-beta filter for velocity.
   * Should be called in a loop (e.g., in opcontrol).
   */
  void update();

  /**
   * @brief Resets motor encoders and IMU rotation.
   */
  void reset_sensors();

  /**
   * @brief Configure the Alpha-Beta filter parameters.
   *
   * @param alpha Position tracking gain (0.0 to 1.0).
   * @param beta Velocity tracking gain (0.0 to 1.0).
   */
  void set_filter_params(double alpha, double beta) {
    this->alpha = alpha;
    this->beta = beta;
  }

private:
  std::vector<std::unique_ptr<pros::v5::Motor>> left_motors;
  std::vector<std::unique_ptr<pros::v5::Motor>> right_motors;
  std::vector<std::unique_ptr<pros::v5::Imu>> imus;
  std::vector<std::unique_ptr<pros::v5::Distance>> distances;

  std::vector<MotorData> left_data;
  std::vector<MotorData> right_data;
  std::vector<double> last_distances;

  // alpha-beta filter parameters
  double alpha = 0.5;
  double beta = 0.1;

  double gear_ratio;
  uint32_t last_time = 0;

  // last known values to handle disconnections
  double last_left_pos = 0;
  double last_right_pos = 0;
  double last_left_vel = 0;
  double last_right_vel = 0;
  double last_heading = 0;
};
