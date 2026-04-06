#include "main.h"
#include "chassis.hpp"
#include <hpipm_common.h>
#include "liblvgl/lvgl.h"
#include "motions/motion.hpp"
#include "motions/sysid.hpp"
#include "mpcc/acados_solver_mpcc.h"
#include "odom/odom.hpp"
#include <vector>

// --- Configuration ---
// Adapt these ports to your physical robot
std::vector<int> LEFT_MOTOR_PORTS = {11, 13};
std::vector<int> RIGHT_MOTOR_PORTS = {-12, -1};
std::vector<int> IMU_PORTS = {3};
std::vector<int> DISTANCE_SENSOR_PORTS = {};
double GEAR_RATIO = 1; // 1:1
pros::v5::MotorGears MOTOR_CARTRIDGE = pros::v5::MotorGears::green;

// --- Global Objects ---
Chassis chassis(LEFT_MOTOR_PORTS, RIGHT_MOTOR_PORTS, IMU_PORTS,
                DISTANCE_SENSOR_PORTS, GEAR_RATIO, 4, 11.5, MOTOR_CARTRIDGE,
                EKFConfig{
                    .Q_v = 0.070366,
                    .Q_w = 0.153613,
                    .R_v_enc = 0.000100,
                    .R_w_enc = 0.000100,
                    .R_w_imu = 0.000662,
                    .slip_tolerance = 9.0,
                });

// Particle Filter setup
std::vector<SensorConfig> sc = {};

PFConfig pfc = {.distance_sensor_stddev = 0.83f,
                .motion_distance_stddev = 0.001f,
                .kinetic_friction = 90.0f,
                .max_drift_speed = 41.0f,
                .motion_heading_stddev = 0.001f,
                .uniform_noise_probability = 1.0f / 78.74016f,
                .random_particle_probability = 0.005f,
                .resampling_threshold = 1.5f};

ParticleFilter pf(500, sc, pfc, 0, 0, 0); // Start at (0,0,0)

MotionManager motions(chassis, pf);

// --- LVGL UI Objects ---
lv_obj_t *status_label = nullptr;
lv_obj_t *odom_x_label = nullptr;
lv_obj_t *odom_y_label = nullptr;
lv_obj_t *odom_theta_label = nullptr;
lv_obj_t *l_vel_label = nullptr;
lv_obj_t *r_vel_label = nullptr;

// --- Task Functions ---
pros::Task *motion_task = nullptr;
lv_timer_t *ui_timer = nullptr;

// LVGL Timer Callback for Thread-Safe Dashboard Updates
void ui_timer_cb(lv_timer_t *timer) {
  Estimate est = pf.estimate();
  char buf[64];
  if (odom_x_label) {
    std::snprintf(buf, sizeof(buf), "X: %.2f", est.x);
    lv_label_set_text(odom_x_label, buf);
  }
  if (odom_y_label) {
    std::snprintf(buf, sizeof(buf), "Y: %.2f", est.y);
    lv_label_set_text(odom_y_label, buf);
  }
  if (odom_theta_label) {
    std::snprintf(buf, sizeof(buf), "Theta: %.2f", est.theta_deg);
    lv_label_set_text(odom_theta_label, buf);
  }

  // Update velocities in top right corner
  if (l_vel_label) {
    std::snprintf(buf, sizeof(buf), "L: %.1f", chassis.get_left_velocity());
    lv_label_set_text(l_vel_label, buf);
  }
  if (r_vel_label) {
    std::snprintf(buf, sizeof(buf), "R: %.1f", chassis.get_right_velocity());
    lv_label_set_text(r_vel_label, buf);
  }
}

void motion_task_fn(void *param) {
  MotionManager *m = static_cast<MotionManager *>(param);
  m->loop();
}

/**
 * Runs initialization code. This occurs as soon as the program is started.
 */
void initialize() {
  chassis.reset_sensors();

  chassis.set_gains({
      .alpha_v = 2.361798,         // m/s²/V
      .alpha_w = 17.800235,        // rad/s²/V
      .lambda_v = 20.000000,       // 1/s
      .lambda_w = 20.000000,       // 1/s
      .max_accel_wall_v = 2.3618,  // m/s²  @ 12 V
      .max_accel_wall_w = 17.8002, // rad/s² @ 12 V
      .kS_v = 1.2074,
      .kS_w = 1.9308,
  });

  chassis.set_stsmc_gains({.Lambda = {1.0, 1.0}, // Lower sensitivity to noise
                           .K1 = {2.0, 2},     // Soften the twist
                           .K2 = {1.5, 1.5},     // Stable integral action
                           .s_eps = 0.5,         // Wider smoothing window
                           .z_max = {12, 12},    // Standard voltage range
                           .Imax_per_side = 5.0});

  // Set up the map for localization
  std::vector<MapSegment> segs = {{-70.5, -70.5, 70.5, -70.5},
                                  {70.5, -70.5, 70.5, 70.5},
                                  {70.5, 70.5, -70.5, 70.5},
                                  {-70.5, 70.5, -70.5, -70.5}};
  pf.set_map(segs);

  // Start the unified Odom + Motion task
  motion_task = new pros::Task(motion_task_fn, &motions, "Motion Task");

  lv_lock();

  ui_timer = lv_timer_create(ui_timer_cb, 50, nullptr);

  // dark background
  lv_obj_set_style_bg_color(lv_screen_active(), lv_color_hex(0x1a1a1a), 0);

  status_label = lv_label_create(lv_screen_active());
  lv_obj_set_style_text_font(status_label, &lv_font_montserrat_20, 0);
  lv_obj_set_style_text_color(status_label, lv_color_hex(0x666666), 0);
  lv_label_set_text(status_label, "886Y Yabadabadoo");
  lv_obj_align(status_label, LV_ALIGN_TOP_MID, 0, 5);

  const lv_color_t purple = lv_color_hex(0xBF94FF);

  odom_x_label = lv_label_create(lv_screen_active());
  lv_obj_set_style_text_font(odom_x_label, &lv_font_montserrat_48, 0);
  lv_obj_set_style_text_color(odom_x_label, purple, 0);
  lv_obj_align(odom_x_label, LV_ALIGN_TOP_LEFT, 20, 40);

  odom_y_label = lv_label_create(lv_screen_active());
  lv_obj_set_style_text_font(odom_y_label, &lv_font_montserrat_48, 0);
  lv_obj_set_style_text_color(odom_y_label, purple, 0);
  lv_obj_align(odom_y_label, LV_ALIGN_TOP_LEFT, 20, 100);

  odom_theta_label = lv_label_create(lv_screen_active());
  lv_obj_set_style_text_font(odom_theta_label, &lv_font_montserrat_48, 0);
  lv_obj_set_style_text_color(odom_theta_label, purple, 0);
  lv_obj_align(odom_theta_label, LV_ALIGN_TOP_LEFT, 20, 160);

  // Velocity indicators (top right corner below team name area)
  l_vel_label = lv_label_create(lv_screen_active());
  lv_obj_set_style_text_font(l_vel_label, &lv_font_montserrat_16, 0);
  lv_obj_set_style_text_color(l_vel_label, lv_palette_main(LV_PALETTE_GREY), 0);
  lv_obj_align(l_vel_label, LV_ALIGN_TOP_RIGHT, -10, 35);

  r_vel_label = lv_label_create(lv_screen_active());
  lv_obj_set_style_text_font(r_vel_label, &lv_font_montserrat_16, 0);
  lv_obj_set_style_text_color(r_vel_label, lv_palette_main(LV_PALETTE_GREY), 0);
  lv_obj_align(r_vel_label, LV_ALIGN_TOP_RIGHT, -10, 55);

  lv_unlock();
}

/**
 * Runs the user autonomous code.
 */
void autonomous() {
  // Outer-Loop Cascade PIDs (output target velocity)
  PIDCoeffs linear_pid = {4.0, 0, 0.5}; // P=4 (10in error -> 40 in/s target)
  PIDCoeffs angular_pid = {4.0, 0,
                           0.5}; // P=4 (30deg error -> 120 deg/s target)

  lv_lock();
  if (status_label) {
    lv_label_set_text(status_label, "Starting Auto Sequence...");
  }
  lv_unlock();

  // 1. Move to a point
  PointTarget target1;
  target1.x = 24.0;
  target1.y = 24.0;
  target1.linear_pid = linear_pid;
  target1.steering_gain = 5;
  target1.max_v = 50.0; // 50 in/s
  target1.slew = 2.0;   // 2 in/s per 10ms (200 in/s^2)
  target1.exit_turn_dist = 7;

  motions.moveTo(target1);
  motions.waitUntilPoint(24, 24, 5);
  motions.cancelMotion();

  // 2. Turn to face a point
  // motions.turnToPoint(0, 0, TurnTarget::Direction::NEAREST, angular_pid,
  // 360.0,
  //                     8.0); // max 360 deg/s, slew 8 deg/s per tick
  // motions.waitUntilAngularVelocity(1.0, false);

  // // 3. Move back to origin
  // PointTarget origin;
  // origin.x = 0;
  // origin.y = 0;
  // origin.linear_pid = linear_pid;
  // origin.steering_gain = 1.5;
  // origin.exit_turn_dist = 4.0;
  // origin.max_v = 50.0;
  // origin.reversed = true;

  // motions.moveTo(origin);
  // motions.waitUntilPoint(0, 0, 1.0);

  // motions.cancelMotion();
}

/**
 * Runs the operator control code.
 */
void opcontrol() {
  pros::Controller master(pros::E_CONTROLLER_MASTER);

  if (motion_task) {
    motion_task->remove();
    delete motion_task;
    motion_task = nullptr;
    pros::delay(100);
  }

  // Kill the dedicated UI timer if you want to stop screen updates in
  // opcontrol if (ui_timer) {
  //   lv_lock();
  //   lv_timer_delete(ui_timer);
  //   ui_timer = nullptr;
  //   lv_unlock();
  // }

  motions.cancelMotion();
  

  // run_linear_sysid(chassis);
  // run_angular_sysid(chassis);

  std::cout << "killed" << std::endl;

  static float prev_heading = chassis.get_rotation() * DEG_TO_RAD;

  while (true) {
    // updated chassis and odom
    chassis.update();

    float curr_heading = chassis.get_rotation() * DEG_TO_RAD;
    auto readings = chassis.get_distance_readings();
    std::vector<float> f_readings(readings.begin(), readings.end());

    pf.update(curr_heading, prev_heading, chassis.get_left_velocity(),
              chassis.get_right_velocity(), f_readings.data(), 0.02);
    prev_heading = curr_heading;

    int left = master.get_analog(ANALOG_LEFT_Y);
    int right = master.get_analog(ANALOG_RIGHT_Y);

    chassis.tank(left * 94.488189, right * 94.488189);

    pros::delay(20);
  }
}