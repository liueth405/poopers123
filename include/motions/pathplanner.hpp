#pragma once


#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>

struct Pose {
  double x;
  double y;
  double theta;
	double curvature = 0.0; // Optional: path curvature at this point
};

// Internal ultra-fast 2D vector (Zero allocations, inlineable)
struct PathVec2 {
  double x{0.0};
  double y{0.0};
  PathVec2() = default;
  PathVec2(double _x, double _y) : x(_x), y(_y) {}
  PathVec2 operator+(const PathVec2& o) const { return {x + o.x, y + o.y}; }
  PathVec2 operator-(const PathVec2& o) const { return {x - o.x, y - o.y}; }
  PathVec2 operator*(double s) const { return {x * s, y * s}; }
  PathVec2 operator/(double s) const { return {x / s, y / s}; }
  PathVec2& operator+=(const PathVec2& o) { x += o.x; y += o.y; return *this; }
  double norm() const { return std::sqrt(x * x + y * y); }
};
inline PathVec2 operator*(double s, const PathVec2& v) { return {v.x * s, v.y * s}; }
 // Waypoint with optional heading flip (for reverse approach/exit)
 struct PathWaypoint {
 double x{0.0};
 double y{0.0};
 bool flip_heading{false}; // If true, tangent direction is reversed 180°
 PathWaypoint() = default;
 PathWaypoint(double _x, double _y, bool _flip = false) : x(_x), y(_y), flip_heading(_flip) {}
 };

class PathPlanner {
private:
  struct BezierSegment {
    PathVec2 p0, p1, p2, p3;
    double length;
    PathVec2 evaluate(double t) const;
    PathVec2 derivative(double t) const;
PathVec2 secondDerivative(double t) const;
  };

  std::vector<BezierSegment> segments_;
  std::vector<double> accumulated_length_;
 std::vector<bool> segment_flip_; // Flip heading for each segment
  double total_length_;

  double computeArcLengthGauss(const BezierSegment &seg, double a, double b) const;
  void buildSegments(const std::vector<PathVec2> &waypoints);
 void buildSegments(const std::vector<PathWaypoint> &waypoints);
  Pose computePose(int segment_idx, double t) const;

public:
  PathPlanner() : total_length_(0.0) {}
  
  // Takes YOUR BlasfeoVec format
  PathPlanner(const std::vector<PathVec2> &waypoints) {
    setWaypoints(waypoints);
  }
  
  // Takes YOUR BlasfeoVec format
  void setWaypoints(const std::vector<PathVec2> &waypoints) {
    if (waypoints.size() < 2) {
      throw std::invalid_argument("Need at least 2 waypoints");
    }
    // Convert BlasfeoVec to ultra-fast Vec2 instantly
    std::vector<PathVec2> fast_waypoints;
    for (PathVec2 wp : waypoints) {
      fast_waypoints.push_back(wp);
    }
    buildSegments(fast_waypoints);
  }

 // Takes PathWaypoint with flip flags
 void setWaypoints(const std::vector<PathWaypoint> &waypoints) {
 if (waypoints.size() < 2) {
 throw std::invalid_argument("Need at least 2 waypoints");
 }
 buildSegments(waypoints);
 }

  double getTotalLength() const { return total_length_; }
double getCurvature(double s) const; // Get curvature at arc length s
double projectFromPosition(double x, double y, double s_hint) const; // Project position to path
  Pose getPose(double s) const;
  Pose getPose(double s, int& last_index) const; // Optimized O(1) lookup
};