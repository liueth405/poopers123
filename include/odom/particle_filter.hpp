#ifndef PARTICLE_FILTER_HPP
#define PARTICLE_FILTER_HPP

#include "approximation_method.hpp"
#include "chassis.h"
#include "exp.h"
#include "normal.h"
#include <algorithm>
#include <arm_neon.h>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <random>
#include <vector>


constexpr float MAP_HALF = 70.5f;
constexpr float MAP_MIN = -MAP_HALF;
constexpr float MAP_MAX = MAP_HALF;
constexpr float PI = 3.14159265358979323846f;
constexpr float PI_2 = PI / 2.0f;
constexpr float PI2 = 2 * PI;
constexpr float RAD_TO_DEG = 180.0f / PI;
constexpr float DEG_TO_RAD = PI / 180.0f;
constexpr float INV_SQRT_2PI = 0.3989422804014327f;
constexpr float MAX_DIST_RECIPROCAL = 1.0f / 78.74016f;
// Use inline constexpr to avoid ODR/linker issues from header-defined globals.
inline constexpr float goldenRatio = 0.61803398875f; // (sqrt(5)-1)/2
inline constexpr float PI_SQRT3 = -1.81379936423f;   // -PI/sqrt(3)

template <typename T> T clamp(T val, T mn, T mx) {
  return std::max(std::min(val, mx), mn);
}

struct PFInput {
  float compass_heading;
  float prev_heading;
  float vl;
  float vr;
  float lateralAccel;
  float *readings;
  float dt;
  float lateralDist;
};

struct Particle {
  float x, y, theta, lateralVel, weight, cos_a, sin_a;
};

struct SensorConfig {
  float dx, dy, angle;
};

class ParticleFilter {
  std::vector<Particle> particles;
  std::mt19937 gen{std::random_device{}()};
  std::ranlux24_base de;
  float cum_distance = 0.0f;
  float alpha_rand = 0.005f;
  float phi = 0.0203177062516f; // should be 0.0203177062516f for real fields
  float stddev = 0.83f;
  float stddevNoise = .45f;
  float stdDrift = 0.01f; // low for drift
  float prand = 1 / 78.74016f;
  float trackwidth = 13.5;
  float f_k = 0.169474562 * 386.2205 * 0.01;

  // Additional parameters to ensure unit consistency and reduce drift:
  float resample_distance_threshold; // threshold (in map units) to
                                     // trigger resampling

  // Sensor cache (local robot frame)
  struct Sensor {
    float dx, dy; // sensor offset: dx forward, dy right
    float angle;  // sensor mounting offset (in radians, compass convention)
    float cos_a, sin_a;
  };

  // wrap angle to prevent extra large numbers
  void anglewrap(float &angle) const {
    if (angle > PI) {
      angle -= PI2;
    } else if (angle < -PI) {
      angle += PI2;
    }
  }

  template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }

  std::vector<Sensor> sensors;

public:
  // The constructor now accepts optional parameters for scaling and noise
  // tuning.
  float cos_h, sin_h;
  float prevV_x = 0;
  ParticleFilter(size_t n, const std::vector<SensorConfig> &sc, float start_x,
                 float start_y, float theta, float start_spread = 2.0f,
                 float noise = 0.25f,
                 float resample_distance_threshold_ = 0.01f) {
    resample_distance_threshold = resample_distance_threshold_;
    particles.resize(n);
    std::normal_distribution<float> dist_x(start_x, start_spread);
    std::normal_distribution<float> dist_y(start_y, start_spread);
    normal_setup();

    for (auto &p : particles) {
      p.x = clamp(dist_x(gen), MAP_MIN, MAP_MAX);
      p.y = clamp(dist_y(gen), MAP_MIN, MAP_MAX);
      p.weight = 1.0f / n;
      p.theta = normal() * 0.05f + theta;
      anglewrap(p.theta);
      p.cos_a = std::cos(p.theta);
      p.sin_a = std::sin(p.theta);
      p.lateralVel = 0.0f;
    }
    // Convert sensor mounting angles from degrees to radians.
    for (const auto &s : sc) {
      float angle = s.angle * DEG_TO_RAD;
      sensors.push_back(
          {s.dx, s.dy, angle, std::cos(PI_2 - angle), std::sin(PI_2 - angle)});
    }
  }

  void reset(double x, double y, float theta) {
    std::normal_distribution<float> dist_x(x, 1);
    std::normal_distribution<float> dist_y(y, 1);

    for (auto &p : particles) {
      p.x = clamp(dist_x(gen), MAP_MIN, MAP_MAX);
      p.y = clamp(dist_y(gen), MAP_MIN, MAP_MAX);
      p.theta = normal() * 0.005f + theta;
      anglewrap(p.theta);
      p.lateralVel = 0.0f;
      p.weight = 1.0f / 2500.0f;
    }
  }

  // precondition: compass_heading: in radians, prev_heading in radians
  void update(float compass_heading, float prev_heading, float vl, float vr,
              float lateralAccel, const float *readings, float dt,
              float lateralDist = 0.0f) {
    float distanceBase = (vr + vl) / 2.0f;
    cum_distance += std::abs(distanceBase);
    for (auto &p : particles) {
      float distance = distanceBase + normal() * stddevNoise;
      const float a_lat = lateralAccel + normal() * stdDrift;

      const float v_lat_new = p.lateralVel + a_lat * dt;
      float lateralDis = 0.5f * (p.lateralVel + v_lat_new) * dt;

      float delta_heading =
          compass_heading - prev_heading +
          normal() * 0.01f; // TODO: you need an actual variable for ts

      // half-step
      float ch = std::cos(0.5f * delta_heading);
      float sh = std::sin(0.5f * delta_heading);

      float c = 2.0f * ch * ch - 1.0f;
      float s = 2.0f * sh * ch;

      float ca_new = p.cos_a * c - p.sin_a * s;
      float sa_new = p.sin_a * c + p.cos_a * s;

      p.cos_a = ca_new;
      p.sin_a = sa_new;

      float cosH = p.cos_a * ch - p.sin_a * sh;
      float sinH = p.sin_a * ch + p.cos_a * sh;

      float deltaY = 0.0f;
      if (std::fabs(delta_heading) > 1e-6f) {
        // deltaY = 2.0f * std::sin(delta_heading * 0.5f) * (distance /
        // delta_heading);
        deltaY = 2.0f * sh * (distance / delta_heading);
      } else {
        deltaY = distance;
      }
      float deltaX = lateralDist != 0 ? lateralDist : lateralDis;

      // Update particle pose in field frame
      p.x = clamp(p.x + deltaX * cosH + deltaY * sinH, MAP_MIN, MAP_MAX);
      p.y = clamp(p.y + -deltaX * sinH + deltaY * cosH, MAP_MIN, MAP_MAX);

      p.theta += delta_heading;
      anglewrap(p.theta);

      p.lateralVel = v_lat_new;
    }

    // Sensor update: adjust particle weights based on sensor readings.
    size_t size_sensor = (sensors.empty() ? 0 : sensors.size() - 1);
    float total_weight = 0.0f;
    for (auto &p : particles) {
      float weight = 1.0f;
      // cos_h = std::cos(p.theta);
      // sin_h = std::sin(p.theta);
      for (size_t i = 0; i < sensors.size(); ++i) {
        const auto &s = sensors[i];
        // Skip invalid/unavailable sensor measurements
        const float meas = readings[i];
        if (!std::isfinite(meas) || meas >= 9998.0f) {
          continue;
        }

        // float x = s.dx * cos_h + s.dy * sin_h;
        // float y = -s.dx * sin_h + s.dy * cos_h;
        float x = s.dx * p.cos_a + s.dy * p.sin_a;
        float y = -s.dx * p.sin_a + s.dy * p.cos_a;

        // TODO: use the angle addition formula to avoid recomputing sin and cos
        // float cos_s = std::cos(PI_2 - (p.theta + s.angle));
        // float sin_s = std::sin(PI_2 - (p.theta + s.angle));
        float cos_s = s.cos_a * p.cos_a + s.sin_a * p.sin_a;
        float sin_s = s.sin_a * p.cos_a - s.cos_a * p.sin_a;

        const float w_i = calc_weight(p.x + x, p.y + y, cos_s, sin_s, meas);
      }
      if (!std::isfinite(weight) || weight <= 0.0f)
        weight = 1e-6f;
      p.weight = weight;
      total_weight += p.weight;
    }
    // Normalize weights with guards
    if (!std::isfinite(total_weight) || total_weight <= 0.0f) {
      const float uniform = 1.0f / std::max<size_t>(particles.size(), 1);
      for (auto &p : particles)
        p.weight = uniform;
    } else {
      const float inv_total = 1.0f / total_weight;
      for (auto &p : particles)
        p.weight *= inv_total;
    }
    // Resample if enough (scaled) distance has been accumulated.
    if (cum_distance < resample_distance_threshold)
      return;
    cum_distance = 0.0f;
    resample();
  }

  void resample() {
    // Systematic resampling
    const size_t numParticles = particles.size();
    if (numParticles > 0) {
      std::vector<float> cdf(numParticles);
      float acc = 0.0f;
      for (size_t i = 0; i < numParticles; ++i) {
        acc += particles[i].weight;
        cdf[i] = acc;
      }
      cdf[numParticles - 1] = 1.0f;

      std::uniform_real_distribution<float> u0dist(
          0.0f, 1.0f / static_cast<float>(numParticles));
      const float u0 = u0dist(gen);
      const float step = 1.0f / static_cast<float>(numParticles);

      std::vector<Particle> oldParticles = particles;
      size_t j = 0;
      for (size_t i = 0; i < numParticles; ++i) {
        const float u = u0 + i * step;
        while (j + 1 < numParticles && u > cdf[j])
          ++j;
        particles[i] = oldParticles[j];
        particles[i].weight = step;
      }
    }
  }

  std::vector<float> estimate() {
    float x = 0.0f, y = 0.0f, theta = 0.0f, lateral = 0.0f, c = 0.0f, s = 0.0f;
    for (auto &p : particles) {
      x += p.x * p.weight;
      y += p.y * p.weight;
      lateral += p.lateralVel * p.weight;

      // s += std::sin(p.theta) * p.weight;
      // c += std::cos(p.theta) * p.weight;
      // TODO: use p.cos_a and p.sin_a to avoid recomputing
      s += p.sin_a * p.weight;
      c += p.cos_a * p.weight;
    }
    x = clamp(x, MAP_MIN, MAP_MAX);
    y = clamp(y, MAP_MIN, MAP_MAX);
    theta = clamp(atan2f(s, c), -PI, PI);
    return std::vector<float>{x, y, theta * RAD_TO_DEG, lateral};
  }

private:
  // Compute the distance from (x0, y0) to the nearest wall along the direction
  // theta. theta is given in math coordinates (0 along positive X, CCW
  // positive).
  float wall_distance(double x0, double y0, float cos_theta, float sin_theta) {
    const double wall_min = MAP_MIN;
    const double wall_max = MAP_MAX;
    double min_dist = std::numeric_limits<double>::infinity();

    // Check intersection with vertical walls.
    if (fabs(cos_theta) > 1e-6) {
      double target_x = (cos_theta < 0) ? wall_min : wall_max;
      double t_x = (target_x - x0) / cos_theta;
      if (t_x >= 0) {
        double y_inter = y0 + t_x * sin_theta;
        if (y_inter >= wall_min && y_inter <= wall_max)
          min_dist = t_x;
      }
    }

    // Check intersection with horizontal walls.
    if (fabs(sin_theta) > 1e-6) {
      double target_y = (sin_theta < 0) ? wall_min : wall_max;
      double t_y = (target_y - y0) / sin_theta;
      if (t_y >= 0) {
        double x_inter = x0 + t_y * cos_theta;
        if (x_inter >= wall_min && x_inter <= wall_max)
          min_dist = std::min(min_dist, t_y);
      }
    }
    return min_dist;
  }

  float calc_weight(float x, float y, float cos_theta, float sin_theta,
                    float measurement) {
    // Compute the expected distance to the wall from the robot's pose and
    // orientation.
    // z from distance sensor
    // d = distance / expected distance
    float distance = wall_distance(x, y, cos_theta, sin_theta);

    // Normalize the difference by the sensor noise (stddev).
    float diff = (measurement - distance) / stddev;

    float32x4_t exp_values = fast_exp_neon(
        (float32x4_t){PI_SQRT3 * (distance - 77.5f) / stddev, -phi * distance,
                      -phi * measurement, -0.5f * diff * diff});
    float weight = 0;
    float alpha_randSub = 1 - alpha_rand;
    float alpha_hit = vgetq_lane_f32(exp_values, 1) * alpha_randSub;

    if (distance <= 77.5f) {
      weight += prand * alpha_rand + alpha_hit * (INV_SQRT_2PI / stddev) *
                                         vgetq_lane_f32(exp_values, 3);
    }

    // p short
    if (measurement < distance) {
      weight += alpha_randSub * phi * vgetq_lane_f32(exp_values, 2);
    }

    // p_max
    if (measurement >= 77.5f) {
      weight += alpha_hit * (1 / (1 + vgetq_lane_f32(exp_values, 0)));
    }

    // Return the final combined weight.
    return weight;
  }
};

/* should be replaced soon with odom.hpp implementation
namespace localization
{
  ParticleFilter pf(1500, std::vector<SensorConfig>({{0, 0, 270}, {0, 0, 0}, {0,
0, 180}}), 0, 0, 0); Eigen::Vector3d state; pros::Distance left(1);
  pros::Distance right(2);
  pros::Distance back(3);
  double prev_left = 0;
  double prev_right = 0;
  double prev_heading = 0;
  void setup(double left, double right, double heading)
  {
    prev_left = left;
    prev_right = right;
    prev_heading = heading;
  }
  void update_robot_state(Chassis chassis)
  {
    float leftDistance = 70;
    float rightDistance = 9999;
    float backDistance = 75;

    float measurements[3] = {leftDistance, rightDistance, backDistance};

    double delta_left = chassis.leftPosition() - prev_left;
    double delta_right = chassis.rightPosition() - prev_right;

    pf.update(chassis.getHeading() * M_PI / 180, prev_heading, delta_left,
delta_right, 0, measurements, 0.01); auto statePf = pf.estimate();

    state << statePf[0], statePf[1], statePf[2];

    prev_left = chassis.leftPosition();
    prev_right = chassis.rightPosition();
    prev_heading = chassis.getHeading() * M_PI / 180;
  }
  Eigen::Vector3d get_robot_state()
  {
    return state;
  }
}*/

#endif // PARTICLE_FILTER_HPP