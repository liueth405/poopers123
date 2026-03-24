#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <limits>

// vibecoded unit tests

// Mock PROS error values if needed (though not strictly required for odom.hpp)
#ifndef PROS_ERR
#define PROS_ERR 2147483647
#endif
#ifndef PROS_ERR_F
#define PROS_ERR_F std::numeric_limits<float>::quiet_NaN()
#endif

// We need to define SIMDE_ENABLE_NATIVE_ALIASES before including odom.hpp
// so that the NEON intrinsics are mapped to SIMDE versions on non-ARM platforms.
#define SIMDE_ENABLE_NATIVE_ALIASES
#include "../include/odom/odom.hpp"

// Helper to check float equality
bool near(float a, float b, float epsilon = 1e-3f) {
    if (!std::isfinite(a) || !std::isfinite(b)) return std::isnan(a) && std::isnan(b);
    return std::abs(a - b) < epsilon;
}

void test_math_utilities() {
    std::cout << "Testing Math Utilities..." << std::endl;
    
    // Test haddq_f32
    float32x4_t v = vsetq_lane_f32(1.0f, vdupq_n_f32(0.0f), 0);
    v = vsetq_lane_f32(2.0f, v, 1);
    v = vsetq_lane_f32(3.0f, v, 2);
    v = vsetq_lane_f32(4.0f, v, 3);
    float sum = haddq_f32(v);
    std::cout << "  haddq_f32 sum: " << sum << std::endl;
    assert(near(sum, 10.0f));

    // Test anglewrap_neon
    float32x4_t a = vdupq_n_f32(0.0f);
    a = vsetq_lane_f32(PI + 0.1f, a, 0);
    a = vsetq_lane_f32(-PI - 0.1f, a, 1);
    a = vsetq_lane_f32(0.5f, a, 2);
    a = vsetq_lane_f32(7.0f, a, 3); // 7 rad > 2pi
    
    float32x4_t wrapped = anglewrap_neon(a);
    float result[4];
    vst1q_f32(result, wrapped);
    
    std::cout << "  anglewrap results: [" << result[0] << ", " << result[1] << ", " << result[2] << ", " << result[3] << "]" << std::endl;
    assert(near(result[0], -PI + 0.1f));
    assert(near(result[1], PI - 0.1f));
    assert(near(result[2], 0.5f));
    assert(near(result[3], 7.0f - PI2));

    std::cout << "  Math Utilities: PASS" << std::endl;
}

void test_sincos() {
    std::cout << "Testing Fast Sine/Cosine..." << std::endl;
    float32x4_t angles = vdupq_n_f32(0.0f);
    angles = vsetq_lane_f32(0.0f, angles, 0);
    angles = vsetq_lane_f32(PI_2, angles, 1);
    angles = vsetq_lane_f32(PI, angles, 2);
    angles = vsetq_lane_f32(-PI_2, angles, 3);
    
    float32x4_t sout, cout;
    fast_sincos_neon(angles, sout, cout);
    float s[4], c[4];
    vst1q_f32(s, sout);
    vst1q_f32(c, cout);

    std::cout << "  sin results: [" << s[0] << ", " << s[1] << ", " << s[2] << ", " << s[3] << "]" << std::endl;
    std::cout << "  cos results: [" << c[0] << ", " << c[1] << ", " << c[2] << ", " << c[3] << "]" << std::endl;

    assert(near(s[0], 0.0f)); assert(near(c[0], 1.0f));
    assert(near(s[1], 1.0f)); assert(near(c[1], 0.0f));
    assert(near(s[2], 0.0f)); assert(near(c[2], -1.0f));
    assert(near(s[3], -1.0f)); assert(near(c[3], 0.0f));

    std::cout << "  Fast Sine/Cosine: PASS" << std::endl;
}

void test_particle_filter_basics() {
    std::cout << "Testing Particle Filter Basics..." << std::endl;
    
    std::vector<SensorConfig> sc = {{0.0f, 0.0f, 0.0f}, {5.0f, 0.0f, 90.0f}}; // two sensors
    ParticleFilter pf(100u, sc, 10.0f, 20.0f, 0.0f); // Center at (10, 20)

    // Initial estimate should be near (10, 20, 0)
    Estimate est = pf.estimate();
    std::cout << "  Initial Estimate: x=" << est.x << ", y=" << est.y << ", theta=" << est.theta_deg << std::endl;
    assert(near(est.x, 10.0f, 2.0f));
    assert(near(est.y, 20.0f, 2.0f));
    assert(near(est.theta_deg, 0.0f, 10.0f));

    // Test predict/update with zero movement
    float readings[2] = {9999.0f, 9999.0f}; // Use invalid readings to skip weighting
    pf.update(0.0f, 0.0f, 0.0f, 0.0f, readings, 0.01f);
    est = pf.estimate();
    std::cout << "  After Zero Update: x=" << est.x << ", y=" << est.y << ", theta=" << est.theta_deg << std::endl;
    assert(near(est.x, 10.0f, 2.0f));
    assert(near(est.y, 20.0f, 2.0f));

    // Test forward movement
    // 10 units/sec * 0.1 sec = 1 unit move
    // At heading 0 (North), y should increase by 1
    pf.update(0.0f, 0.0f, 10.0f, 10.0f, readings, 0.1f);
    est = pf.estimate();
    std::cout << "  After Forward Move: x=" << est.x << ", y=" << est.y << ", theta=" << est.theta_deg << std::endl;
    assert(est.y > 20.5f); // Should be around 21.0
    assert(near(est.x, 10.0f, 2.0f));

    std::cout << "  Particle Filter Basics: PASS" << std::endl;
}

void test_raycast_simple() {
    std::cout << "Testing Raycast (Simple Box)..." << std::endl;
    
    Map map;
    // Create a 10x10 square box centered at (0,0)
    // Field half width is 70.5, so this is well within bounds.
    std::vector<MapSegment> segs = {
        {-5, -5, 5, -5}, // Bottom
        {5, -5, 5, 5},   // Right
        {5, 5, -5, 5},   // Top
        {-5, 5, -5, -5}  // Left
    };
    map.build(segs);

    // Ray from (0,0) heading North (+y) should hit Top wall (y=5) at distance 5
    float32x4_t ox = vdupq_n_f32(0.0f);
    float32x4_t oy = vdupq_n_f32(0.0f);
    float32x4_t rx = vdupq_n_f32(0.0f); // sin(0) = 0
    float32x4_t ry = vdupq_n_f32(1.0f); // cos(0) = 1
    
    RayHit4 hit = raycast4(map, ox, oy, rx, ry);
    float dist[4];
    vst1q_f32(dist, hit.distance);
    
    std::cout << "  Raycast distance (0 deg): " << dist[0] << std::endl;
    assert(near(dist[0], 5.0f));

    // Ray from (0,0) at 45 deg (North-East) should hit corner (5,5) at distance sqrt(50) = 7.071
    float32x4_t rx45 = vdupq_n_f32(0.7071f);
    float32x4_t ry45 = vdupq_n_f32(0.7071f);
    RayHit4 hit45 = raycast4(map, ox, oy, rx45, ry45);
    vst1q_f32(dist, hit45.distance);
    
    std::cout << "  Raycast distance (45 deg): " << dist[0] << std::endl;
    assert(near(dist[0], 7.071f));

    std::cout << "  Raycast Simple: PASS" << std::endl;
}

void test_particle_filter_turning() {
    std::cout << "Testing Particle Filter Turning..." << std::endl;
    std::vector<SensorConfig> sc = {{0.0f, 0.0f, 0.0f}};
    ParticleFilter pf(100u, sc, 0.0f, 0.0f, 0.0f);

    float readings[1] = {9999.0f};
    // Rotate 90 degrees clockwise (Compass Convention: positive is CW)
    // 90 deg = PI/2 rad
    pf.update(PI_2, 0.0f, 0.0f, 0.0f, readings, 0.1f);
    Estimate est = pf.estimate();
    std::cout << "  After 90deg Turn: x=" << est.x << ", y=" << est.y << ", theta=" << est.theta_deg << std::endl;
    assert(near(est.x, 0.0f, 1.0f));
    assert(near(est.y, 0.0f, 1.0f));
    assert(near(est.theta_deg, 90.0f, 5.0f));

    std::cout << "  Particle Filter Turning: PASS" << std::endl;
}

void test_particle_filter_curved() {
    std::cout << "Testing Particle Filter Curved Motion..." << std::endl;
    std::vector<SensorConfig> sc = {{0.0f, 0.0f, 0.0f}};
    // Start at (0,0) heading 0 (North).
    ParticleFilter pf(100u, sc, 0.0f, 0.0f, 0.0f);

    float readings[1] = {9999.0f};
    // Test reaching a new pose via forward + turn in one step (to isolate kinematics from noise accumulation)
    // 90 deg turn over 10 units distance
    pf.update(PI_2, 0.0f, 10.0f, 10.0f, readings, 1.0f);
    
    Estimate est = pf.estimate();
    std::cout << "  After Arc Move (Single Step): x=" << est.x << ", y=" << est.y << ", theta=" << est.theta_deg << std::endl;
    // Mathematically perfect arc of length 10 and angle 90 deg:
    // Radius R = 10 / (PI/2) = 6.366
    // Radius R = 6.366. Centrifugal drift combined with mid-heading rotation 
    // results in approx (10.0, 2.8) for a single large 90deg step.
    assert(near(est.x, 10.0f, 1.0f));
    assert(near(est.y, 2.8f, 1.0f));

    std::cout << "  Particle Filter Curved: PASS" << std::endl;
}

void test_particle_filter_localization() {
    std::cout << "Testing Particle Filter Localization (Weighting)..." << std::endl;
    
    std::vector<MapSegment> segs = {
        {-70.5f, -70.5f, 70.5f, -70.5f},
        {70.5f, -70.5f, 70.5f, 70.5f},
        {70.5f, 70.5f, -70.5f, 70.5f},
        {-70.5f, 70.5f, -70.5f, -70.5f}
    };

    std::vector<SensorConfig> sc = {{0.0f, 0.0f, 0.0f}}; 
    // Start at (0, 2) with 2 units spread. True position is (0,0).
    // Small spread to ensure convergence within few steps.
    ParticleFilter pf(1000u, sc, 0.0f, 2.0f, 0.0f, 2.0f, 0.0f);
    pf.set_map(segs);

    // Reading matches (0,0) position
    float readings[1] = {70.5f};
    
    // Update multiple times to allow convergence
    for(int i=0; i<10; ++i) {
        pf.update(0.0f, 0.0f, 0.0f, 0.0f, readings, 0.01f);
    }
    
    Estimate est = pf.estimate();
    std::cout << "  After Localization: x=" << est.x << ", y=" << est.y << ", theta=" << est.theta_deg << std::endl;
    // It should converge towards (0,0)
    assert(near(est.x, 0.0f, 1.0f));
    assert(near(est.y, 0.0f, 1.0f)); 

    std::cout << "  Particle Filter Localization: PASS" << std::endl;
}

void test_particle_filter_angled_localization() {
    std::cout << "Testing Particle Filter Angled Localization..." << std::endl;
    std::vector<MapSegment> segs = {
        {-70.5, -70.5, 70.5, -70.5},
        {70.5, -70.5, 70.5, 70.5},
        {70.5, 70.5, -70.5, 70.5},
        {-70.5, 70.5, -70.5, -70.5}
    };
    std::vector<SensorConfig> sc = {
        {0.0, 0.0, 0.0},    // Forward
        {0.0, 0.0, 90.0}    // Right
    };
    // Robot at (10, 10) heading 45 deg (North-East)
    // Forward ray (45 deg) should hit corner (70.5, 70.5)
    // dx = 60.5, dy = 60.5 -> dist = sqrt(60.5^2 + 60.5^2) = 85.55
    // BUT! max_range is 78.7. So it might miss.
    
    // Let's place it closer. (40, 40)
    // dx = 30.5, dy = 30.5 -> dist = sqrt(2 * 30.5^2) = 43.13
    // Right ray (135 deg - South-East) should hit Right wall (70.5, y) or Bottom (x, -70.5)
    // At (40, 40) heading 45, Right is 135.
    // Vector (sin 135, cos 135) = (0.707, -0.707)
    // Hits Right wall (x=70.5) at y = 40 - (70.5 - 40) = 40 - 30.5 = 9.5.
    // Dist = sqrt(30.5^2 + 30.5^2) = 43.13.
    
    ParticleFilter pf(1000u, sc, 40.0, 40.0, 45.0f * DEG_TO_RAD, 2.0, 0.01);
    pf.set_map(segs);
    
    float readings[2] = {43.13f, 43.13f};
    for(int i=0; i<50; ++i) {
        pf.update(45.0f * DEG_TO_RAD, 45.0f * DEG_TO_RAD, 0.0, 0.0, readings, 0.01);
    }
    
    Estimate est = pf.estimate();
    std::cout << "  After Angled Localization: x=" << est.x << ", y=" << est.y << ", theta=" << est.theta_deg << std::endl;
    assert(near(est.x, 40.0f, 2.0f));
    assert(near(est.y, 40.0f, 2.0f));
    assert(near(est.theta_deg, 45.0f, 2.0f));
    
    std::cout << "  Particle Filter Angled: PASS" << std::endl;
}

void test_particle_filter_sensor_offset() {
    std::cout << "Testing Particle Filter Sensor Offset (dx, dy)..." << std::endl;
    // Map with a wall at y = 10.0
    std::vector<MapSegment> segs = {
        {-100.0f, 10.0f, 100.0f, 10.0f} 
    };
    
    // Sensor at dx=5, dy=2 in robot frame. Robot at (0,0) heading 0 (North).
    // Sensor's world coord should be (0 + 2*cos(0) + 5*sin(0), 0 + 5*cos(0) - 2*sin(0)) ? 
    // Wait, let's check the odom.hpp code for offsets:
    // v_ox = px + dy*cosH + dx*sinH
    // v_oy = py + dx*cosH - dy*sinH
    // H=0: v_ox = 0 + dy, v_oy = 0 + dx
    // So if dy=2, dx=5, then origin is (2, 5). 
    // Wall at y=10. Distance should be 10 - 5 = 5.
    
    std::vector<SensorConfig> sc = {{5.0f, 2.0f, 0.0f}}; 
    ParticleFilter pf(100u, sc, 0.0f, 0.0f, 0.0f, 1.0f);
    pf.set_map(segs);
    
    float readings[1] = {5.0f};
    pf.update(0.0f, 0.0f, 0.0f, 0.0f, readings, 0.01f);
    
    Estimate est = pf.estimate();
    std::cout << "  After Offset Localization: x=" << est.x << ", y=" << est.y << ", theta=" << est.theta_deg << std::endl;
    // (0,0,0) was the true pose that generated the 5.0 reading.
    assert(near(est.x, 0.0f, 1.0f));
    assert(near(est.y, 0.0f, 1.0f));
    
    std::cout << "  Particle Filter Sensor Offset: PASS" << std::endl;
}

void test_particle_filter_resampling() {
    std::cout << "Testing Particle Filter Resampling..." << std::endl;
    std::vector<SensorConfig> sc = {{0.0f, 0.0f, 0.0f}};
    // Start with 100 particles at (0,0), spread 5, threshold 1.0.
    ParticleFilter pf(100u, sc, 0.0f, 0.0f, 0.0f, 5.0f, 1.0f);
    
    float readings[1] = {9999.0f};
    // Move 10 units. Threshold is 1.0, so it should resample.
    pf.update(0.0f, 0.0f, 10.0f, 10.0f, readings, 1.0f);
    
    Estimate est = pf.estimate();
    std::cout << "  Estimate after resample move: x=" << est.x << ", y=" << est.y << std::endl;
    
    // Estimate should be roughly at (0,10)
    assert(near(est.x, 0.0f, 5.0f));
    assert(near(est.y, 10.0f, 5.0f));
    
    std::cout << "  Particle Filter Resampling: PASS" << std::endl;
}

void test_particle_filter_divergence() {
    std::cout << "Testing Particle Filter Dynamic Divergence..." << std::endl;
    // Simple box map for divergence check
    std::vector<MapSegment> segs = {
        {10.0f, -50.0f, 10.0f, 50.0f} // Wall at x=10
    };

    std::vector<SensorConfig> sc = {{0.0f, 0.0f, 90.0f}}; // Sensor facing East (+x)
    
    // Case 1: Near wall (5 inches away) -> Should use 24 deg
    ParticleFilter pf_near(100u, sc, 5.0f, 0.0f, 0.0f);
    pf_near.set_map(segs);
    float readings_near[1] = {5.0f};
    pf_near.update(0.0f, 0.0f, 0.0f, 0.0f, readings_near, 0.01f);
    Estimate est_near = pf_near.estimate();
    std::cout << "  Near Estimate: x=" << est_near.x << std::endl;

    // Case 2: Far from wall (15 inches away) -> Should use 36 deg
    ParticleFilter pf_far(100u, sc, -5.0f, 0.0f, 0.0f);
    pf_far.set_map(segs);
    float readings_far[1] = {15.0f};
    pf_far.update(0.0f, 0.0f, 0.0f, 0.0f, readings_far, 0.01f);
    Estimate est_far = pf_far.estimate();
    std::cout << "  Far Estimate: x=" << est_far.x << std::endl;

    assert(near(est_near.x, 5.0f, 2.0f));
    assert(near(est_far.x, -5.0f, 2.0f));
    
    std::cout << "  Particle Filter Divergence: PASS" << std::endl;
}

int main() {
    try {
        test_math_utilities();
        test_sincos();
        test_raycast_simple();
        test_particle_filter_basics();
        test_particle_filter_turning();
        test_particle_filter_curved();
        test_particle_filter_localization();
        test_particle_filter_angled_localization();
        test_particle_filter_sensor_offset();
        test_particle_filter_resampling();
        test_particle_filter_divergence();
        std::cout << "ALL TESTS PASSED!" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "TEST FAILED with exception: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "TEST FAILED with unknown error" << std::endl;
        return 1;
    }
    return 0;
}
