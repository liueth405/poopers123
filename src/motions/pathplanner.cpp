#include "motions/pathplanner.hpp"
#include <limits>

PathVec2 PathPlanner::BezierSegment::evaluate(double t) const {
  double u = 1.0 - t, tt = t * t, uu = u * u;
  double uuu = uu * u, ttt = tt * t;
  PathVec2 p = p0 * uuu;
  p += p1 * (3.0 * uu * t);
  p += p2 * (3.0 * u * tt);
  p += p3 * ttt;
  return p;
}

PathVec2 PathPlanner::BezierSegment::derivative(double t) const {
  double u = 1.0 - t;
  PathVec2 dp = (p1 - p0) * (3.0 * u * u);
  dp += (p2 - p1) * (6.0 * u * t);
  dp += (p3 - p2) * (3.0 * t * t);
  return dp;
}

Pose PathPlanner::getPose(double s, int &last_index) const {
  if (segments_.empty())
    return {0.0, 0.0, 0.0};
  if (s <= 0.0) {
    last_index = 0;
    return computePose(0, 0.0);
  }
  if (s >= total_length_) {
    last_index = segments_.size() - 1;
    return computePose(last_index, 1.0);
  }

  // Sequential search starting from last_index (O(1) in steady state)
  int segment_idx = std::clamp(last_index, 0, (int)segments_.size() - 1);

  if (s < accumulated_length_[segment_idx]) {
    // Search backward
    while (segment_idx > 0 && s < accumulated_length_[segment_idx]) {
      segment_idx--;
    }
  } else {
    // Search forward
    while (segment_idx < (int)segments_.size() - 1 &&
           s >= accumulated_length_[segment_idx + 1]) {
      segment_idx++;
    }
  }

  last_index = segment_idx;
  double s_target = s - accumulated_length_[segment_idx];
  const auto &seg = segments_[segment_idx];
  double t = s_target / seg.length;

  // Refine t using simple Newton step (Arc length integration)
  for (int i = 0; i < 4; ++i) {
    double err = computeArcLengthGauss(seg, 0.0, t) - s_target;
    double dt = seg.derivative(t).norm();
    if (std::abs(dt) < 1e-6)
      break;
    t -= err / dt;
    t = std::clamp(t, 0.0, 1.0);
  }

  return computePose(segment_idx, t);
}

Pose PathPlanner::getPose(double s) const {
  int dummy = 0;
  return getPose(s, dummy);
}

double PathPlanner::computeArcLengthGauss(const BezierSegment &seg, double a,
                                          double b) const {
  static const double x[] = {0.0, 0.5384693101056831, -0.5384693101056831,
                             0.9061798459386640, -0.9061798459386640};
  static const double w[] = {0.5688888888888889, 0.4786286704993665,
                             0.4786286704993665, 0.2369268850561891,
                             0.2369268850561891};
  double sum = 0.0, hd = 0.5 * (b - a), hs = 0.5 * (a + b);
  for (int i = 0; i < 5; ++i)
    sum += w[i] * seg.derivative(hd * x[i] + hs).norm();
  return hd * sum;
}

void PathPlanner::buildSegments(const std::vector<PathVec2> &waypoints) {
  int n = waypoints.size() - 1;
  segments_.resize(n);
  segment_flip_.clear(); // No flip flags for PathVec2 version
  if (n < 1)
    return;

  std::vector<double> L(n);
  for (int i = 0; i < n; ++i)
    L[i] = std::max((waypoints[i + 1] - waypoints[i]).norm(), 1e-6);

  if (n == 1) {
    segments_[0].p0 = waypoints[0];
    segments_[0].p1 =
        waypoints[0] + (waypoints[1] - waypoints[0]) * (1.0 / 3.0);
    segments_[0].p2 =
        waypoints[0] + (waypoints[1] - waypoints[0]) * (2.0 / 3.0);
    segments_[0].p3 = waypoints[1];
  } else {
    // Fast O(N) Thomas Algorithm (No Matrix Inversions, No Heap Overloads)
    std::vector<double> cp(n + 1);
    std::vector<PathVec2> dp(n + 1), V(n + 1);

    cp[0] = 1.0 / 2.0;
    dp[0] = (waypoints[1] - waypoints[0]) * (3.0 / L[0]) / 2.0;

    for (int i = 1; i < n; ++i) {
      double denom = 2.0 * (L[i - 1] + L[i]) - L[i] * cp[i - 1];
      cp[i] = L[i - 1] / denom;
      PathVec2 d_i =
          (waypoints[i + 1] - waypoints[i]) * (3.0 * L[i - 1] / L[i]) +
          (waypoints[i] - waypoints[i - 1]) * (3.0 * L[i] / L[i - 1]);
      dp[i] = (d_i - dp[i - 1] * L[i]) / denom;
    }

    double denom_n = 2.0 - cp[n - 1];
    dp[n] = ((waypoints[n] - waypoints[n - 1]) * (3.0 / L[n - 1]) - dp[n - 1]) /
            denom_n;

    V[n] = dp[n];
    for (int i = n - 1; i >= 0; --i)
      V[i] = dp[i] - V[i + 1] * cp[i];

    for (int i = 0; i < n; ++i) {
      segments_[i].p0 = waypoints[i];
      segments_[i].p1 = waypoints[i] + V[i] * (L[i] / 3.0);
      segments_[i].p2 = waypoints[i + 1] - V[i + 1] * (L[i] / 3.0);
      segments_[i].p3 = waypoints[i + 1];
    }
  }

  accumulated_length_.assign(1, 0.0);
  for (int i = 0; i < n; ++i) {
    segments_[i].length = computeArcLengthGauss(segments_[i], 0.0, 1.0);
    accumulated_length_.push_back(accumulated_length_.back() +
                                  segments_[i].length);
  }
  total_length_ = accumulated_length_.back();
}

Pose PathPlanner::computePose(int segment_idx, double t) const {
  PathVec2 p = segments_[segment_idx].evaluate(t);
  PathVec2 dp = segments_[segment_idx].derivative(t);
  double theta = std::atan2(dp.x, dp.y);

  // Apply 180 degree flip if this segment has the flip flag set (for reverse
  // motion)
  if (!segment_flip_.empty() && segment_flip_[segment_idx]) {
    theta += M_PI;
    while (theta > M_PI)
      theta -= 2.0 * M_PI;
    while (theta < -M_PI)
      theta += 2.0 * M_PI;
  }

  PathVec2 ddp = segments_[segment_idx].secondDerivative(t);
	double num = std::abs(dp.x * ddp.y - dp.y * ddp.x);
	double denom = std::pow(dp.x * dp.x + dp.y * dp.y, 1.5);
	double curvature = (denom < 1e-10) ? 0.0 : num / denom;
	return {p.x, p.y, theta, curvature};
}

// Overload that handles flip flags
void PathPlanner::buildSegments(const std::vector<PathWaypoint> &waypoints) {
  int n = waypoints.size() - 1;
  segments_.resize(n);
  segment_flip_.resize(n);
  if (n < 1)
    return;

  std::vector<double> L(n);
  for (int i = 0; i < n; ++i) {
    L[i] = std::max(std::sqrt((waypoints[i + 1].x - waypoints[i].x) *
                                  (waypoints[i + 1].x - waypoints[i].x) +
                              (waypoints[i + 1].y - waypoints[i].y) *
                                  (waypoints[i + 1].y - waypoints[i].y)),
                    1e-6);
    // Store flip flag at the END of each segment
    segment_flip_[i] = waypoints[i + 1].flip_heading;
  }

  // Same spline construction but using PathWaypoint
  std::vector<double> cp(n + 1);
  std::vector<double> dp_x(n + 1), dp_y(n + 1);
  std::vector<double> V_x(n + 1), V_y(n + 1);

  cp[0] = 1.0 / 2.0;
  dp_x[0] = (waypoints[1].x - waypoints[0].x) * (3.0 / L[0]) / 2.0;
  dp_y[0] = (waypoints[1].y - waypoints[0].y) * (3.0 / L[0]) / 2.0;

  for (int i = 1; i < n; ++i) {
    double denom = 2.0 * (L[i - 1] + L[i]) - L[i] * cp[i - 1];
    cp[i] = L[i - 1] / denom;
    double d_x =
        (waypoints[i + 1].x - waypoints[i].x) * (3.0 * L[i - 1] / L[i]) +
        (waypoints[i].x - waypoints[i - 1].x) * (3.0 * L[i] / L[i - 1]);
    double d_y =
        (waypoints[i + 1].y - waypoints[i].y) * (3.0 * L[i - 1] / L[i]) +
        (waypoints[i].y - waypoints[i - 1].y) * (3.0 * L[i] / L[i - 1]);
    dp_x[i] = (d_x - dp_x[i - 1] * L[i]) / denom;
    dp_y[i] = (d_y - dp_y[i - 1] * L[i]) / denom;
  }

  double denom_n = 2.0 - cp[n - 1];
  dp_x[n] =
      ((waypoints[n].x - waypoints[n - 1].x) * (3.0 / L[n - 1]) - dp_x[n - 1]) /
      denom_n;
  dp_y[n] =
      ((waypoints[n].y - waypoints[n - 1].y) * (3.0 / L[n - 1]) - dp_y[n - 1]) /
      denom_n;

  V_x[n] = dp_x[n];
  V_y[n] = dp_y[n];
  for (int i = n - 1; i >= 0; --i) {
    V_x[i] = dp_x[i] - V_x[i + 1] * cp[i];
    V_y[i] = dp_y[i] - V_y[i + 1] * cp[i];
  }

  for (int i = 0; i < n; ++i) {
    segments_[i].p0 = {waypoints[i].x, waypoints[i].y};
    segments_[i].p1 = {waypoints[i].x + V_x[i] * (L[i] / 3.0),
                       waypoints[i].y + V_y[i] * (L[i] / 3.0)};
    segments_[i].p2 = {waypoints[i + 1].x - V_x[i + 1] * (L[i] / 3.0),
                       waypoints[i + 1].y - V_y[i + 1] * (L[i] / 3.0)};
    segments_[i].p3 = {waypoints[i + 1].x, waypoints[i + 1].y};
  }

  accumulated_length_.assign(1, 0.0);
  for (int i = 0; i < n; ++i) {
    segments_[i].length = computeArcLengthGauss(segments_[i], 0.0, 1.0);
    accumulated_length_.push_back(accumulated_length_.back() +
                                  segments_[i].length);
  }
  total_length_ = accumulated_length_.back();
}

PathVec2 PathPlanner::BezierSegment::secondDerivative(double t) const {
    double u = 1.0 - t;
    PathVec2 ddp = (p2 - p1 - (p1 - p0)) * (6.0 * u);
    ddp += (p3 - p2 - (p2 - p1)) * (6.0 * t);
    return ddp;
}

double PathPlanner::getCurvature(double s) const {
    if (segments_.empty()) return 0.0;
    
    s = std::clamp(s, 0.0, total_length_);
    
    // Find segment
    int segment_idx = 0;
    for (size_t i = 0; i < segments_.size(); ++i) {
        if (s < accumulated_length_[i + 1]) {
            segment_idx = i;
            break;
        }
    }
    
    double s_local = s - accumulated_length_[segment_idx];
    const auto &seg = segments_[segment_idx];
    
    // Newton's method to find t
    double t = s_local / seg.length;
    for (int iter = 0; iter < 4; ++iter) {
        double err = computeArcLengthGauss(seg, 0.0, t) - s_local;
        double dt = seg.derivative(t).norm();
        if (std::abs(dt) < 1e-6) break;
        t -= err / dt;
        t = std::clamp(t, 0.0, 1.0);
    }
    
    // Curvature: κ = |x'y'' - y'x''| / (x'² + y'²)^(3/2)
    PathVec2 d1 = seg.derivative(t);
    PathVec2 d2 = seg.secondDerivative(t);
    
    double numerator = std::abs(d1.x * d2.y - d1.y * d2.x);
    double denominator = std::pow(d1.x * d1.x + d1.y * d1.y, 1.5);
    
    if (denominator < 1e-10) return 0.0;
    
    return numerator / denominator;
}

double PathPlanner::projectFromPosition(double x, double y, double s_hint) const {
    if (segments_.empty()) return 0.0;
    
    // Start search from hint
    double s_best = std::clamp(s_hint, 0.0, total_length_);
    double dist_best = std::numeric_limits<double>::max();
    
    // Search nearby segment
    int segment_idx = 0;
    for (size_t i = 0; i < segments_.size(); ++i) {
        if (s_best < accumulated_length_[i + 1]) {
            segment_idx = i;
            break;
        }
    }
    
    // Newton iteration to minimize distance to path
    double s = s_best;
    for (int iter = 0; iter < 10; ++iter) {
        Pose pose = getPose(s);
        double dx = x - pose.x;
        double dy = y - pose.y;
        
        // Distance squared
        double dist_sq = dx * dx + dy * dy;
        if (dist_sq < dist_best) {
            dist_best = dist_sq;
            s_best = s;
        }
        
        // Gradient: d(dist)/ds ≈ -tangent · (point - pose)
        // tangent = (sin(theta), cos(theta)) for compass heading
        double tangent_x = std::sin(pose.theta);
        double tangent_y = std::cos(pose.theta);
        double grad = -(tangent_x * dx + tangent_y * dy);
        
        // Hessian approximation (second derivative of distance)
        double hess = 1.0; // Simplified: assume unit curvature
        
        double ds = -grad / hess;
        s += ds * 0.5; // Damped update
        
        // Clamp to valid range
        s = std::clamp(s, 0.0, total_length_);
        
        if (std::abs(ds) < 1e-4) break;
    }
    
    return s_best;
}
