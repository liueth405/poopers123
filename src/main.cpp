#include "main.h"
#include "odom/odom.hpp"
#include "motions/motion.hpp"
#include <vector>

// --- Configuration ---
// Adaptive these ports to your physical robot
std::vector<int8_t> LEFT_MOTOR_PORTS = {1, 2, 3};
std::vector<int8_t> RIGHT_MOTOR_PORTS = {-4, -5, -6};
std::vector<int8_t> IMU_PORTS = {7};
std::vector<int8_t> DISTANCE_SENSOR_PORTS = {8, 9};
double GEAR_RATIO = 0.6; // 36:60
pros::v5::MotorGears MOTOR_CARTRIDGE = pros::v5::MotorGears::blue;

// --- Global Objects ---
Chassis chassis(LEFT_MOTOR_PORTS, RIGHT_MOTOR_PORTS, IMU_PORTS, DISTANCE_SENSOR_PORTS, GEAR_RATIO, MOTOR_CARTRIDGE);

// Particle Filter setup
std::vector<SensorConfig> sc = {
    {0.0, 0.0, 0.0}, // Forward distance sensor
    {0.0, 0.0, 90.0} // Right distance sensor
};

PFConfig pfc = {
    .distance_sensor_stddev = 0.83f,
    .motion_distance_stddev = 0.45f,
    .lateral_viscous_friction = 0.5f,
    .motion_heading_stddev = 0.01f,
    .uniform_noise_probability = 1.0f / 78.74016f,
    .random_particle_probability = 0.005f,
    .resampling_threshold = 0.01f
};

ParticleFilter pf(500, sc, pfc, 0, 0, 0); // Start at (0,0,0)

MotionManager motions(chassis, pf);

// --- Task Function ---
pros::Task* motion_task = nullptr;
void motion_task_fn(void* param) {
    MotionManager* m = static_cast<MotionManager*>(param);
    m->loop();
}

/**
 * Runs initialization code. This occurs as soon as the program is started.
 */
void initialize() {
    pros::lcd::initialize();
    
    // Set up the map for localization
    std::vector<MapSegment> segs = {
        {-70.5, -70.5, 70.5, -70.5},
        {70.5, -70.5, 70.5, 70.5},
        {70.5, 70.5, -70.5, 70.5},
        {-70.5, 70.5, -70.5, -70.5}
    };
    pf.set_map(segs);

    // Start the unified Odom + Motion task
    motion_task = new pros::Task(motion_task_fn, &motions, "Motion Task");
    
    pros::lcd::set_text(1, "Chassis & Motions Initialized");
}

/**
 * Runs the user autonomous code.
 */
void autonomous() {
    // Example sequence using the new Motion Wrapper
    PIDCoeffs linear_pid = {100, 0, 10};
    PIDCoeffs angular_pid = {50, 0, 5};

    pros::lcd::set_text(2, "Starting Auto Sequence...");

    // 1. Move to a point
    PointTarget target1;
    target1.x = 24.0;
    target1.y = 24.0;
    target1.linear_pid = linear_pid;
    target1.angular_pid = angular_pid;
    target1.max_v = 8000;
    target1.slew = 400;
    
    motions.moveTo(target1);
    motions.waitUntilPoint(24, 24, 2.0);
    
    // 2. Turn to face a point
    motions.turnToPoint(0, 0, TurnTarget::Direction::NEAREST, angular_pid, 6000, 300);
    motions.waitUntilAngle(180, 2.0); // Facing (0,0) from (24,24) is roughly 225 deg? 
    // Actually turnToPoint calculates it, so we can just waitUntilDone
    while(!motions.isDone()) pros::delay(20);

    // 3. Move back to origin
    PointTarget origin;
    origin.x = 0;
    origin.y = 0;
    origin.linear_pid = linear_pid;
    origin.angular_pid = angular_pid;
    origin.max_v = 10000;
    
    motions.moveTo(origin);
    motions.waitUntilPoint(0, 0, 1.0);

    motions.cancelMotion();


}

/**
 * Runs the operator control code.
 */
void opcontrol() {
    pros::Controller master(pros::E_CONTROLLER_MASTER);
    
    // Kill the background motion/odom task to prevent interference
    if (motion_task) {
        motion_task->remove();
        delete motion_task;
        motion_task = nullptr;
    }
    
    motions.cancelMotion();

    static float prev_heading = 0;

    while (true) {
        // Since we killed the background task, we must update sensors manually 
        // if we still want Odom telemetry in opcontrol.
        chassis.update();
        
        float curr_heading = chassis.get_rotation() * DEG_TO_RAD;
        auto readings = chassis.get_distance_readings();
        std::vector<float> f_readings(readings.begin(), readings.end());
        pf.update(curr_heading, prev_heading, chassis.get_left_velocity(), chassis.get_right_velocity(), f_readings.data(), 0.02);
        prev_heading = curr_heading;

        // Tank drive control
        int left = master.get_analog(ANALOG_LEFT_Y);
        int right = master.get_analog(ANALOG_RIGHT_Y);
        
        // Convert -127 to 127 range to -12000 to 12000 range
        chassis.tank(left * 94, right * 94);

        // Display current pose from Odom
        Estimate est = pf.estimate();
        pros::lcd::print(4, "X: %.2f  Y: %.2f", est.x, est.y);
        pros::lcd::print(5, "Theta: %.2f", est.theta_deg);

        pros::delay(20);
    }
}