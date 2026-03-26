#ifndef POOPERS_ODOM_HPP_FINAL
#define POOPERS_ODOM_HPP_FINAL

#include <algorithm>
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
#else
#define SIMDE_ENABLE_NATIVE_ALIASES
#include "simde/arm/neon.h"
#endif
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <random>
#include <vector>

#include "exp.h"
#include "normal.h"

constexpr float MAP_HALF = 70.5f;
constexpr float MAP_MIN = -MAP_HALF;
constexpr float MAP_MAX = MAP_HALF;
constexpr float PI = 3.14159265358979323846f;
constexpr float PI_2 = PI / 2.0f;
constexpr float PI2 = 2.0f * PI;
constexpr float RAD_TO_DEG = 180.0f / PI;
constexpr float DEG_TO_RAD = PI / 180.0f;
constexpr float INV_SQRT_2PI = 0.3989422804014327f;
inline constexpr float PI_SQRT3 = -1.81379936423f; // -π / √3

struct SensorConfig {
  float dx, dy, angle;
};

struct Sensor {
  float dx, dy;
  float angle;
  float cos_a, sin_a;
  // pre-computed ray orientations for 24 deg [0] and 36 deg [1] divergence
  // x = cos(angle + phi), y = sin(angle + phi)
  float ray_sin[2][8];
  float ray_cos[2][8];
  float ray_w[2][8];
};

struct MapSegment {
  float x1, y1, x2, y2;
};

struct Map {
  size_t num_segments = 0;
  size_t num_aligned = 0;
  float *seg_x1 = nullptr;
  float *seg_y1 = nullptr;
  float *seg_dx = nullptr;
  float *seg_dy = nullptr;
  float *seg_inv_len = nullptr;
  float *mem = nullptr;

  // spatial partitioning (grid)
  static constexpr int GRID_RES = 16;
  uint16_t *grid_heads = nullptr;   // offset into grid_indices
  uint8_t *grid_counts = nullptr;   // count of segments in cell
  uint16_t *grid_indices = nullptr; // flattened index list
  float grid_min_x, grid_min_y;
  float grid_cell_size;
  float grid_inv_cell_size;

  Map() = default;

  void build(const std::vector<MapSegment> &segs) {
    num_segments = segs.size();
    num_aligned = (num_segments + 3u) & ~3u;

    if (mem)
      free(mem);
    if (grid_heads)
      free(grid_heads);
    if (grid_counts)
      free(grid_counts);
    if (grid_indices)
      free(grid_indices);

    size_t total_bytes = 5 * num_aligned * sizeof(float) + 64;
    mem = (float *)malloc(total_bytes);
    uintptr_t addr = reinterpret_cast<uintptr_t>(mem);
    float *aligned_mem =
        reinterpret_cast<float *>((addr + 63u) & ~uintptr_t(63u));

    seg_x1 = aligned_mem + 0 * num_aligned;
    seg_y1 = aligned_mem + 1 * num_aligned;
    seg_dx = aligned_mem + 2 * num_aligned;
    seg_dy = aligned_mem + 3 * num_aligned;
    seg_inv_len = aligned_mem + 4 * num_aligned;

    for (size_t i = 0; i < num_segments; ++i) {
      const MapSegment &s = segs[i];
      float dx = s.x2 - s.x1;
      float dy = s.y2 - s.y1;
      seg_x1[i] = s.x1;
      seg_y1[i] = s.y1;
      seg_dx[i] = dx;
      seg_dy[i] = dy;
      float len = std::sqrt(dx * dx + dy * dy);
      seg_inv_len[i] = (len > 1e-8f) ? 1.0f / len : 0.0f;
    }
    for (size_t i = num_segments; i < num_aligned; ++i)
      seg_x1[i] = seg_y1[i] = seg_dx[i] = seg_dy[i] = seg_inv_len[i] = 0.0f;

    // build grid
    grid_min_x = MAP_MIN;
    grid_min_y = MAP_MIN;
    grid_cell_size = (MAP_MAX - MAP_MIN) / (float)GRID_RES;
    grid_inv_cell_size = 1.0f / grid_cell_size;

    grid_heads = (uint16_t *)calloc(GRID_RES * GRID_RES, sizeof(uint16_t));
    grid_counts = (uint8_t *)calloc(GRID_RES * GRID_RES, sizeof(uint8_t));

    // pass 1: count segments per cell
    for (size_t i = 0; i < num_segments; ++i) {
      float xmin = std::min(seg_x1[i], seg_x1[i] + seg_dx[i]);
      float xmax = std::max(seg_x1[i], seg_x1[i] + seg_dx[i]);
      float ymin = std::min(seg_y1[i], seg_y1[i] + seg_dy[i]);
      float ymax = std::max(seg_y1[i], seg_y1[i] + seg_dy[i]);

      int ix1 = std::max(0, std::min(GRID_RES - 1, (int)((xmin - grid_min_x) *
                                                         grid_inv_cell_size)));
      int ix2 = std::max(0, std::min(GRID_RES - 1, (int)((xmax - grid_min_x) *
                                                         grid_inv_cell_size)));
      int iy1 = std::max(0, std::min(GRID_RES - 1, (int)((ymin - grid_min_y) *
                                                         grid_inv_cell_size)));
      int iy2 = std::max(0, std::min(GRID_RES - 1, (int)((ymax - grid_min_y) *
                                                         grid_inv_cell_size)));

      for (int iy = iy1; iy <= iy2; ++iy) {
        for (int ix = ix1; ix <= ix2; ++ix) {
          int idx = iy * GRID_RES + ix;
          if (grid_counts[idx] < 255)
            grid_counts[idx]++;
        }
      }
    }

    // prefix sum to find heads
    size_t total_indices = 0;
    for (int i = 0; i < GRID_RES * GRID_RES; ++i) {
      grid_heads[i] = (uint16_t)total_indices;
      total_indices += grid_counts[i];
      grid_counts[i] = 0; // reset for second pass
    }

    grid_indices = (uint16_t *)malloc(total_indices * sizeof(uint16_t));

    // pass 2: store indices
    for (size_t i = 0; i < num_segments; ++i) {
      float xmin = std::min(seg_x1[i], seg_x1[i] + seg_dx[i]);
      float xmax = std::max(seg_x1[i], seg_x1[i] + seg_dx[i]);
      float ymin = std::min(seg_y1[i], seg_y1[i] + seg_dy[i]);
      float ymax = std::max(seg_y1[i], seg_y1[i] + seg_dy[i]);

      int ix1 = std::max(0, std::min(GRID_RES - 1, (int)((xmin - grid_min_x) *
                                                         grid_inv_cell_size)));
      int ix2 = std::max(0, std::min(GRID_RES - 1, (int)((xmax - grid_min_x) *
                                                         grid_inv_cell_size)));
      int iy1 = std::max(0, std::min(GRID_RES - 1, (int)((ymin - grid_min_y) *
                                                         grid_inv_cell_size)));
      int iy2 = std::max(0, std::min(GRID_RES - 1, (int)((ymax - grid_min_y) *
                                                         grid_inv_cell_size)));

      for (int iy = iy1; iy <= iy2; ++iy) {
        for (int ix = ix1; ix <= ix2; ++ix) {
          int idx = iy * GRID_RES + ix;
          int pos = grid_heads[idx] + grid_counts[idx];
          grid_indices[pos] = (uint16_t)i;
          grid_counts[idx]++;
        }
      }
    }
  }

  ~Map() {
    if (mem)
      free(mem);
    if (grid_heads)
      free(grid_heads);
    if (grid_counts)
      free(grid_counts);
    if (grid_indices)
      free(grid_indices);
  }
  Map(const Map &) = delete;
  Map &operator=(const Map &) = delete;
};

struct RayHit4 {
  float32x4_t distance;
  float32x4_t normal_x;
  float32x4_t normal_y;
};

struct Estimate {
  float x, y, theta_deg, lateral;
};

template <typename T> static inline T clamp(T val, T mn, T mx) {
  return std::max(std::min(val, mx), mn);
}

// VIBE CODED:
static inline float32x4_t fast_recip(float32x4_t x) {
  float32x4_t r = vrecpeq_f32(x);
  r = vmulq_f32(r, vrecpsq_f32(x, r));
  return r;
}

static inline float haddq_f32(float32x4_t v) {
  float32x2_t s = vadd_f32(vget_low_f32(v), vget_high_f32(v));
  s = vpadd_f32(s, s);
  return vget_lane_f32(s, 0);
}

static inline float32x4_t anglewrap_neon(float32x4_t a) {
  const float32x4_t vpi = vdupq_n_f32(PI);
  const float32x4_t v2pi = vdupq_n_f32(PI2);
  const float32x4_t vnpi = vnegq_f32(vpi);
  uint32x4_t hi = vcgtq_f32(a, vpi);
  uint32x4_t lo = vcltq_f32(a, vnpi);
  a = vsubq_f32(
      a, vreinterpretq_f32_u32(vandq_u32(hi, vreinterpretq_u32_f32(v2pi))));
  a = vaddq_f32(
      a, vreinterpretq_f32_u32(vandq_u32(lo, vreinterpretq_u32_f32(v2pi))));
  return a;
}

// Cephes-style vectorized sin/cos for 4 floats at once
static inline void fast_sincos_neon(float32x4_t x, float32x4_t &out_sin,
                                    float32x4_t &out_cos) {
  const float32x4_t v4op = vdupq_n_f32(1.27323954473516f);
  const float32x4_t vdp1 = vdupq_n_f32(-0.78515625f);
  const float32x4_t vdp2 = vdupq_n_f32(-2.4187564849853515625e-4f);
  const float32x4_t vdp3 = vdupq_n_f32(-3.77489497744594108e-8f);
  const float32x4_t vsp0 = vdupq_n_f32(-1.9515295891e-4f);
  const float32x4_t vsp1 = vdupq_n_f32(8.3321608736e-3f);
  const float32x4_t vsp2 = vdupq_n_f32(-1.6666654611e-1f);
  const float32x4_t vcp0 = vdupq_n_f32(2.443315711809948e-5f);
  const float32x4_t vcp1 = vdupq_n_f32(-1.388731625493765e-3f);
  const float32x4_t vcp2 = vdupq_n_f32(4.166664568298827e-2f);
  const float32x4_t vhalf = vdupq_n_f32(0.5f);
  const float32x4_t vone = vdupq_n_f32(1.0f);

  uint32x4_t sign_sin =
      vandq_u32(vreinterpretq_u32_f32(x), vdupq_n_u32(0x80000000u));
  x = vabsq_f32(x);

  float32x4_t y = vmulq_f32(x, v4op);
  uint32x4_t j = vcvtq_u32_f32(y);
  j = vandq_u32(vaddq_u32(j, vdupq_n_u32(1)), vdupq_n_u32(~1u));
  y = vcvtq_f32_u32(j);

  sign_sin = veorq_u32(sign_sin, vshlq_n_u32(vandq_u32(j, vdupq_n_u32(4)), 29));
  uint32x4_t sign_cos =
      vshlq_n_u32(vbicq_u32(vdupq_n_u32(4), vsubq_u32(j, vdupq_n_u32(2))), 29);

  uint32x4_t poly_mask =
      vceqq_u32(vandq_u32(j, vdupq_n_u32(2)), vdupq_n_u32(0));

  x = vmlaq_f32(x, y, vdp1);
  x = vmlaq_f32(x, y, vdp2);
  x = vmlaq_f32(x, y, vdp3);
  float32x4_t z = vmulq_f32(x, x);

  float32x4_t yc = vmlaq_f32(vcp1, vcp0, z);
  yc = vmlaq_f32(vcp2, yc, z);
  yc = vmulq_f32(yc, vmulq_f32(z, z));
  yc = vmlsq_f32(yc, z, vhalf);
  yc = vaddq_f32(yc, vone);

  float32x4_t ys = vmlaq_f32(vsp1, vsp0, z);
  ys = vmlaq_f32(vsp2, ys, z);
  ys = vmulq_f32(ys, z);
  ys = vmlaq_f32(x, ys, x);

  float32x4_t s_val = vbslq_f32(poly_mask, ys, yc);
  float32x4_t c_val = vbslq_f32(poly_mask, yc, ys);

  out_sin =
      vreinterpretq_f32_u32(veorq_u32(vreinterpretq_u32_f32(s_val), sign_sin));
  out_cos =
      vreinterpretq_f32_u32(veorq_u32(vreinterpretq_u32_f32(c_val), sign_cos));
}

// ── 4-Particle Vectorized Raycast
// ─────────────────────────────────────────────
static inline RayHit4 raycast4(const Map &map, float32x4_t ox, float32x4_t oy,
                               float32x4_t rx, float32x4_t ry) {
  const float32x4_t v_zero = vdupq_n_f32(0.0f);
  const float32x4_t v_one = vdupq_n_f32(1.0f);
  const float32x4_t v_eps = vdupq_n_f32(1e-8f);
  const float32x4_t v_inf = vdupq_n_f32(std::numeric_limits<float>::infinity());
  const float max_range = 70.5f;

  float32x4_t best_t = v_inf;
  float32x4_t best_nx = v_zero;
  float32x4_t best_ny = v_zero;

  // [OPTIMIZATION: Vectorized Ray AABB Generation]
  // Instead of sequentially calculating limits using scalar float arrays,
  // we calculate the maximum end-point vectors exactly in NEON mathematics.
  // vminq/vmaxq give us the axis-aligned bounding box (AABB) bounds for all 4
  // particles simultaneously. This bounding box is then mapped back into
  // scalars to cull the grid cells that the 4 clustered rays definitely cannot
  // touch.
  const float32x4_t v_max_range = vdupq_n_f32(max_range);
  const float32x4_t ex = vmlaq_f32(ox, rx, v_max_range);
  const float32x4_t ey = vmlaq_f32(oy, ry, v_max_range);

  float xmina[4], xmaxa[4], ymina[4], ymaxa[4];
  vst1q_f32(xmina, vminq_f32(ox, ex));
  vst1q_f32(xmaxa, vmaxq_f32(ox, ex));
  vst1q_f32(ymina, vminq_f32(oy, ey));
  vst1q_f32(ymaxa, vmaxq_f32(oy, ey));

  float xmin = std::min({xmina[0], xmina[1], xmina[2], xmina[3]});
  float xmax = std::max({xmaxa[0], xmaxa[1], xmaxa[2], xmaxa[3]});
  float ymin = std::min({ymina[0], ymina[1], ymina[2], ymina[3]});
  float ymax = std::max({ymaxa[0], ymaxa[1], ymaxa[2], ymaxa[3]});

  int ix1 = clamp((int)((xmin - map.grid_min_x) * map.grid_inv_cell_size), 0,
                  Map::GRID_RES - 1);
  int ix2 = clamp((int)((xmax - map.grid_min_x) * map.grid_inv_cell_size), 0,
                  Map::GRID_RES - 1);
  int iy1 = clamp((int)((ymin - map.grid_min_y) * map.grid_inv_cell_size), 0,
                  Map::GRID_RES - 1);
  int iy2 = clamp((int)((ymax - map.grid_min_y) * map.grid_inv_cell_size), 0,
                  Map::GRID_RES - 1);

  uint32_t visited[32] = {0};

  for (int iy = iy1; iy <= iy2; ++iy) {
    for (int ix = ix1; ix <= ix2; ++ix) {
      int idx = iy * Map::GRID_RES + ix;
      int count = map.grid_counts[idx];
      if (!count)
        continue;

      int head = map.grid_heads[idx];
      for (int k = 0; k < count; ++k) {
        uint16_t s_idx = map.grid_indices[head + k];
        if (visited[s_idx >> 5] & (1u << (s_idx & 31)))
          continue;
        visited[s_idx >> 5] |= (1u << (s_idx & 31));

        const float32x4_t sx = vdupq_n_f32(map.seg_x1[s_idx]);
        const float32x4_t sy = vdupq_n_f32(map.seg_y1[s_idx]);
        const float32x4_t sdx = vdupq_n_f32(map.seg_dx[s_idx]);
        const float32x4_t sdy = vdupq_n_f32(map.seg_dy[s_idx]);
        const float32x4_t sinv = vdupq_n_f32(map.seg_inv_len[s_idx]);

        const float32x4_t qpx = vsubq_f32(sx, ox);
        const float32x4_t qpy = vsubq_f32(sy, oy);
        const float32x4_t det = vmlsq_f32(vmulq_f32(rx, sdy), ry, sdx);
        const float32x4_t abs_det = vabsq_f32(det);
        const uint32x4_t det_ok = vcgtq_f32(abs_det, v_eps);
        const float32x4_t safe_det = vbslq_f32(det_ok, det, v_one);
        const float32x4_t inv_det = fast_recip(safe_det);

        const float32x4_t t =
            vmulq_f32(vmlsq_f32(vmulq_f32(qpx, sdy), qpy, sdx), inv_det);
        const float32x4_t s =
            vmulq_f32(vmlsq_f32(vmulq_f32(qpx, ry), qpy, rx), inv_det);

        const uint32x4_t valid = vandq_u32(
            det_ok,
            vandq_u32(vcgeq_f32(t, v_zero),
                      vandq_u32(vcgeq_f32(s, v_zero), vcleq_f32(s, v_one))));

        const uint32x4_t better = vandq_u32(valid, vcltq_f32(t, best_t));
        best_t = vbslq_f32(better, t, best_t);

        // Update normals for better hits
        float32x4_t nx = vnegq_f32(vmulq_f32(sdy, sinv));
        float32x4_t ny = vmulq_f32(sdx, sinv);
        float32x4_t dot_nr = vaddq_f32(vmulq_f32(nx, rx), vmulq_f32(ny, ry));
        uint32x4_t flip = vcgtq_f32(dot_nr, v_zero);
        nx = vbslq_f32(flip, vnegq_f32(nx), nx);
        ny = vbslq_f32(flip, vnegq_f32(ny), ny);
        best_nx = vbslq_f32(better, nx, best_nx);
        best_ny = vbslq_f32(better, ny, best_ny);
      }
    }
  }
  return {best_t, best_nx, best_ny};
}

// ── Particle Filter
// ───────────────────────────────────────────────────────────
class ParticleFilter {
  size_t num_particles;
  size_t num_aligned;

  void *raw_mem = nullptr;
  float *mem_block = nullptr;

  float *px, *py, *pt, *plv, *pw, *pca, *psa;
  float *noise_d, *noise_h;
  float *old_px, *old_py, *old_pt, *old_plv, *old_pw, *old_pca, *old_psa;
  float *cdf_buf;

  Map map;
  std::vector<Sensor> sensors;
  std::mt19937 gen{static_cast<uint32_t>(std::random_device{}())};

  float accumulated_distance = 0.0f;
  float random_particle_probability = 0.005f;
  float phi = 0.0203177062516f;

  // divergence parameters (10 degrees total)
  float div_angle = 10.0f * (PI / 180.0f);
  float div_step = div_angle / 7.0f;
  float sigma_sq = (div_angle * 0.5f) * (div_angle * 0.5f);


  // tune these !!! except for max_range_threshold
  float distance_sensor_stddev = 0.83f;
  float motion_distance_stddev = 0.45f;
  float lateral_viscous_friction = 0.5f;
  float motion_heading_stddev = 0.01f;
  float uniform_noise_probability = 1.0f / 78.74016f;
  float resampling_threshold;

  float max_range = 78.74016f;
  float max_range_threshold = 77.5f;

  static void anglewrap(float &angle) {
    angle = std::fmod(angle + PI, PI2);
    if (angle < 0.0f)
      angle += PI2;
    angle -= PI;
  }

public:
  ParticleFilter(size_t n, const std::vector<SensorConfig> &sc, float start_x,
                 float start_y, float theta, float start_spread = 2.0f,
                 float resampling_threshold_ = 0.01f)
      : num_particles(n), resampling_threshold(resampling_threshold_) {
    num_aligned = (n + 3u) & ~3u;
    size_t total_floats = 17 * num_aligned;
    size_t total_bytes = total_floats * sizeof(float) + 64;
    raw_mem = malloc(total_bytes);
    uintptr_t addr = reinterpret_cast<uintptr_t>(raw_mem);
    mem_block = reinterpret_cast<float *>((addr + 63u) & ~uintptr_t(63u));
    std::memset(mem_block, 0, total_floats * sizeof(float));

    px = mem_block + 0 * num_aligned;
    py = mem_block + 1 * num_aligned;
    pt = mem_block + 2 * num_aligned;
    plv = mem_block + 3 * num_aligned;
    pw = mem_block + 4 * num_aligned;
    pca = mem_block + 5 * num_aligned;
    psa = mem_block + 6 * num_aligned;
    noise_d = mem_block + 7 * num_aligned;
    noise_h = mem_block + 8 * num_aligned;
    old_px = mem_block + 9 * num_aligned;
    old_py = mem_block + 10 * num_aligned;
    old_pt = mem_block + 11 * num_aligned;
    old_plv = mem_block + 12 * num_aligned;
    old_pw = mem_block + 13 * num_aligned;
    old_pca = mem_block + 14 * num_aligned;
    old_psa = mem_block + 15 * num_aligned;
    cdf_buf = mem_block + 16 * num_aligned;

    normal_setup();

    float w0 = 1.0f / (float)n;
    for (size_t i = 0; i < n; ++i) {
      px[i] = clamp(start_x + (float)normal() * start_spread, MAP_MIN, MAP_MAX);
      py[i] = clamp(start_y + (float)normal() * start_spread, MAP_MIN, MAP_MAX);
      float t = theta + (float)normal() * 0.05f;
      anglewrap(t);
      pt[i] = t;
      plv[i] = 0.0f;
      pw[i] = w0;
      pca[i] = std::cos(t);
      psa[i] = std::sin(t);
    }

    for (const auto &s : sc) {
      float a = s.angle * DEG_TO_RAD;
      Sensor ns;
      ns.dx = s.dx;
      ns.dy = s.dy;
      ns.angle = a;
      ns.cos_a = std::cos(a);
      ns.sin_a = std::sin(a);

      float divs[2] = {24.0f * (PI / 180.0f), 36.0f * (PI / 180.0f)};
      for (int m = 0; m < 2; ++m) {
        float d_ang = divs[m];
        float d_step = d_ang / 7.0f;
        float s_sq = (d_ang * 0.5f) * (d_ang * 0.5f);
        for (int r = 0; r < 8; ++r) {
          float phi = -(d_ang * 0.5f) + (float)r * d_step;
          float total_a = a + phi;
          ns.ray_sin[m][r] = std::sin(total_a);
          ns.ray_cos[m][r] = std::cos(total_a);
          ns.ray_w[m][r] = std::exp(-2.0f * (phi * phi) / s_sq);
        }
      }
      sensors.push_back(ns);
    }
    refresh_trig();
  }

  ~ParticleFilter() { free(raw_mem); }
  ParticleFilter(const ParticleFilter &) = delete;
  ParticleFilter &operator=(const ParticleFilter &) = delete;

  void reset(float x, float y, float theta) {
    float w0 = 1.0f / (float)num_particles;
    for (size_t i = 0; i < num_particles; ++i) {
      px[i] = clamp(x + (float)normal() * 1.0f, MAP_MIN, MAP_MAX);
      py[i] = clamp(y + (float)normal() * 1.0f, MAP_MIN, MAP_MAX);
      float t = theta + (float)normal() * 0.005f;
      anglewrap(t);
      pt[i] = t;
      plv[i] = 0.0f;
      pw[i] = w0;
      pca[i] = std::cos(t);
      psa[i] = std::sin(t);
    }
    accumulated_distance = 0.0f;
    refresh_trig();
  }

  void set_map(const std::vector<MapSegment> &segs) { map.build(segs); }

  void update(float compass_heading, float prev_heading, float vl, float vr,
              const float *readings, float dt) {
    predict(compass_heading, prev_heading, vl, vr, dt);

    float total_weight = 0.0f;

    // [OPTIMIZATION: Transcendental Hoisting]
    // The standard deviation and exponential short-hit probabilities only
    // depend on the sensor's static measurement values (vmeas), NOT the
    // particle state. By precalculating `std::exp` perfectly in standard scalar
    // mathematics here, we save the NEON processing pipeline from evaluating
    // `fast_exp_neon` thousands of times inside the inner particle evaluation
    // loop.
    float32x4_t v_inv_sigma = vdupq_n_f32(1.0f / distance_sensor_stddev);
    float32x4_t v_exp_sho_arr[16];
    for (size_t j = 0; j < sensors.size() && j < 16; ++j) {
      v_exp_sho_arr[j] = vdupq_n_f32(std::exp(-phi * readings[j]));
    }

    for (size_t i = 0; i < num_aligned; i += 4) {
      float32x4_t v_weight = vdupq_n_f32(1.0f);
      float32x4_t v_px = vld1q_f32(&px[i]);
      float32x4_t v_py = vld1q_f32(&py[i]);
      float32x4_t v_pca = vld1q_f32(&pca[i]);
      float32x4_t v_psa = vld1q_f32(&psa[i]);

      for (size_t j = 0; j < sensors.size(); ++j) {
        const auto &s = sensors[j];
        const float meas = readings[j];
        if (!std::isfinite(meas) || meas >= 9998.0f)
          continue;

        // compass convention sensor world position (+y 0, CW+)
        // x is east, y is north. s.dx = right, s.dy =fForward
        // E = E_p + s.dx * cos(pt) + s.dy * sin(pt)
        // N = N_p + s.dy * cos(pt) - s.dx * sin(pt)
        float32x4_t v_ox =
            vaddq_f32(v_px, vmlaq_f32(vmulq_f32(vdupq_n_f32(s.dx), v_pca),
                                      vdupq_n_f32(s.dy), v_psa));
        float32x4_t v_oy =
            vaddq_f32(v_py, vmlsq_f32(vmulq_f32(vdupq_n_f32(s.dy), v_pca),
                                      vdupq_n_f32(s.dx), v_psa));

        float32x4_t sw =
            calc_weight4(v_ox, v_oy, v_pca, v_psa, j, vdupq_n_f32(meas),
                         v_inv_sigma, v_exp_sho_arr[j]);
        v_weight = vmulq_f32(v_weight, sw);
      }

      if (i + 4 > num_particles) {
        float mask_arr[4] = {0, 0, 0, 0};
        for (size_t k = 0; k < (num_particles - i); ++k)
          mask_arr[k] = 1.0f;
        v_weight = vmulq_f32(v_weight, vld1q_f32(mask_arr));
      }

      vst1q_f32(&pw[i], v_weight);
      total_weight += haddq_f32(v_weight);
    }

    if (!std::isfinite(total_weight) || total_weight <= 0.0f) {
      float uw = 1.0f / (float)num_particles;
      for (size_t i = 0; i < num_particles; ++i)
        pw[i] = uw;
    } else {
      float32x4_t v_inv = vdupq_n_f32(1.0f / total_weight);
      for (size_t i = 0; i < num_aligned; i += 4) {
        float32x4_t w = vld1q_f32(&pw[i]);
        vst1q_f32(&pw[i], vmulq_f32(w, v_inv));
      }
    }

    if (accumulated_distance < resampling_threshold)
      return;
    accumulated_distance = 0.0f;
    resample();
  }

  Estimate estimate() const {
    float32x4_t sx = vdupq_n_f32(0.0f), sy = vdupq_n_f32(0.0f),
                slv = vdupq_n_f32(0.0f), sc = vdupq_n_f32(0.0f),
                ss = vdupq_n_f32(0.0f);
    for (size_t i = 0; i < num_aligned; i += 4) {
      float32x4_t w = vld1q_f32(&pw[i]);
      sx = vmlaq_f32(sx, vld1q_f32(&px[i]), w);
      sy = vmlaq_f32(sy, vld1q_f32(&py[i]), w);
      slv = vmlaq_f32(slv, vld1q_f32(&plv[i]), w);
      sc = vmlaq_f32(sc, vld1q_f32(&pca[i]), w);
      ss = vmlaq_f32(ss, vld1q_f32(&psa[i]), w);
    }
    return {clamp(haddq_f32(sx), MAP_MIN, MAP_MAX),
            clamp(haddq_f32(sy), MAP_MIN, MAP_MAX),
            atan2f(haddq_f32(ss), haddq_f32(sc)) * RAD_TO_DEG, haddq_f32(slv)};
  }

private:
  void predict(float compass_heading, float prev_heading, float vl, float vr,
               float dt) {
    const float distBase = (vr + vl) * 0.5f;
    const float common_dh = compass_heading - prev_heading;
    accumulated_distance += std::abs(distBase);

    for (size_t i = 0; i < num_particles; ++i) {
      noise_d[i] = (float)normal();
      noise_h[i] = (float)normal();
    }
    for (size_t i = num_particles; i < num_aligned; ++i)
      noise_d[i] = noise_h[i] = 0.0f;

    const float32x4_t v_db = vdupq_n_f32(distBase),
                      v_ms = vdupq_n_f32(motion_distance_stddev);

    const float32x4_t v_cdh = vdupq_n_f32(common_dh),
                      v_hs = vdupq_n_f32(motion_heading_stddev);

    const float32x4_t v_dt = vdupq_n_f32(dt), v_half = vdupq_n_f32(0.5f),
                      v_two = vdupq_n_f32(2.0f), v_one = vdupq_n_f32(1.0f);

    const float32x4_t v_eps = vdupq_n_f32(1e-6f), v_mmin = vdupq_n_f32(MAP_MIN),
                      v_mmax = vdupq_n_f32(MAP_MAX), v_zero = vdupq_n_f32(0.0f);

    for (size_t i = 0; i < num_aligned; i += 4) {
      // get current x, y, theta, lateral velocity, cos(theta_particle), and
      // sin(theta_particle)
      float32x4_t x = vld1q_f32(&px[i]), y = vld1q_f32(&py[i]),
                  t = vld1q_f32(&pt[i]), lv = vld1q_f32(&plv[i]),
                  ca = vld1q_f32(&pca[i]), sa = vld1q_f32(&psa[i]);
      // reopen the noise for distance and heading
      float32x4_t nd = vld1q_f32(&noise_d[i]), nh = vld1q_f32(&noise_h[i]);

      // add the noise to the distance and heading
      float32x4_t dist = vmlaq_f32(v_db, nd, v_ms),
                  dh = vmlaq_f32(v_cdh, nh, v_hs);

      // half-step or delta_theta/2
      float32x4_t half_dh = vmulq_f32(v_half, dh), sh, ch;
      fast_sincos_neon(half_dh, sh, ch);

      // cos and sin double angle identities or delta_theta
      float32x4_t cf = vsubq_f32(vmulq_f32(v_two, vmulq_f32(ch, ch)), v_one),
                  sf = vmulq_f32(v_two, vmulq_f32(sh, ch));

      // update the cos and sin of the particle's heading through cos(theta +
      // delta_theta) = cos(theta)cos(delta_theta) - sin(theta)sin(delta_theta)
      float32x4_t ca_new = vmlsq_f32(vmulq_f32(ca, cf), sa, sf),
                  sa_new = vmlaq_f32(vmulq_f32(sa, cf), ca, sf);

      // update the cos and sin for chord approximation(odom) through cos(theta
      // + delta_theta) = cos(theta)cos(delta_theta/2) -
      // update the cos and sin for chord approximation(odom) through cos(theta
      // + delta_theta) = cos(theta)cos(delta_theta/2) -
      // sin(theta)sin(delta_theta/2)
      float32x4_t cosH = vmlsq_f32(vmulq_f32(ca, ch), sa, sh),
                  sinH = vmlaq_f32(vmulq_f32(sa, ch), ca, sh);

      // handle the case where delta_theta is close to 0
      float32x4_t abs_dh = vabsq_f32(dh);
      // delta_heading > epsilon, deltaY = 2 * sin(delta_theta / 2) /
      // delta_theta
      uint32x4_t small = vcleq_f32(abs_dh, v_eps);
      float32x4_t safe_dh = vbslq_f32(small, v_one, dh),
                  sinc = vbslq_f32(
                      small, v_one,
                      vmulq_f32(vmulq_f32(v_two, sh), fast_recip(safe_dh)));

      // displacement in mid-heading frame (exact chord approximation):
      // in this frame, the lateral component is zero and the forward component
      // is the chord length.
      float32x4_t chord = vmulq_f32(dist, sinc);
      float32x4_t deltaX_kinematic = v_zero;
      float32x4_t deltaY = chord;

      // v_fwd = dist / dt
      float32x4_t v_fwd = vmulq_f32(dist, fast_recip(v_dt));

      // vy_next = (vy - w * dt * v) / (1.0 + (fabs(w) + alpha_p) * dt)
      // w * dt is equivalent to dh (which carries noise).
      // fabs(w) * dt is fabs(dh). alpha_p is lateral_viscous_friction.
      float32x4_t num = vsubq_f32(lv, vmulq_f32(v_fwd, dh));
      float32x4_t den = vaddq_f32(v_one, vaddq_f32(vabsq_f32(dh), vdupq_n_f32(lateral_viscous_friction * dt)));
      float32x4_t lv_raw = vmulq_f32(num, fast_recip(den));
      
      // lv_new = (lv * cf) + (v_fwd * sf) - (lv * k_fric * dt)
      // float32x4_t lv_rot = vaddq_f32(vmulq_f32(lv, cf), vmulq_f32(v_fwd, sf));
      //float32x4_t lv_fric =
      //    vmulq_f32(lv, vdupq_n_f32(lateral_viscous_friction * dt));
      // float32x4_t lv_raw = vsubq_f32(lv_rot, lv_fric);

      // deadband: if lateral velocity < 0.2 inches/s, snap to zero.
      const float32x4_t v_threshold = vdupq_n_f32(0.2f);
      uint32x4_t small_lv = vcltq_f32(vabsq_f32(lv_raw), v_threshold);
      float32x4_t lv_new = vbslq_f32(small_lv, v_zero, lv_raw);

      // add existing slip (lv) to kinematics. d = (v + v0) / 2 * t
      float32x4_t deltaX =
         vaddq_f32(deltaX_kinematic,
                   vmulq_f32(vmulq_f32(v_half, vaddq_f32(lv, lv_new)), v_dt));


      // update the particle's position in field frame (compass convention: +y
      // north, +x east, 0=north, cw+)
      // x' = x + deltaX * cosH + deltaY * sinH
      // y' = y + deltaY * cosH - deltaX * sinH
      float32x4_t nx = vmaxq_f32(
          v_mmin, vminq_f32(v_mmax, vmlaq_f32(vmlaq_f32(x, deltaX, cosH),
                                              deltaY, sinH)));
      float32x4_t ny = vmaxq_f32(
          v_mmin, vminq_f32(v_mmax, vmlsq_f32(vmlaq_f32(y, deltaY, cosH),
                                              deltaX, sinH)));
      float32x4_t nt = anglewrap_neon(vaddq_f32(t, dh));

      // store the new values
      vst1q_f32(&px[i], nx);
      vst1q_f32(&py[i], ny);
      vst1q_f32(&pt[i], nt);
      vst1q_f32(&plv[i], lv_new);
      vst1q_f32(&pca[i], ca_new);
      vst1q_f32(&psa[i], sa_new);
    }
  }

  void resample() {
    if (num_particles == 0)
      return;
    float acc = 0.0f;
    for (size_t i = 0; i < num_particles; ++i) {
      acc += pw[i];
      cdf_buf[i] = acc;
    }
    cdf_buf[num_particles - 1] = 1.0f;
    std::memcpy(old_px, px, num_particles * sizeof(float));
    std::memcpy(old_py, py, num_particles * sizeof(float));
    std::memcpy(old_pt, pt, num_particles * sizeof(float));
    std::memcpy(old_plv, plv, num_particles * sizeof(float));
    std::memcpy(old_pca, pca, num_particles * sizeof(float));
    std::memcpy(old_psa, psa, num_particles * sizeof(float));
    std::uniform_real_distribution<float> u0dist(0.0f,
                                                 1.0f / (float)num_particles);
    const float u0 = u0dist(gen), step = 1.0f / (float)num_particles;
    size_t j = 0;
    for (size_t i = 0; i < num_particles; ++i) {
      const float u = u0 + i * step;
      while (j + 1 < num_particles && u > cdf_buf[j])
        ++j;
      px[i] = old_px[j];
      py[i] = old_py[j];
      pt[i] = old_pt[j];
      plv[i] = old_plv[j];
      pca[i] = old_pca[j];
      psa[i] = old_psa[j];
      pw[i] = step;
    }
    refresh_trig();
  }

  // calc_weight4 assesses the probability of physical sensor measurements
  // matching the simulated map raycasts across a 4-particle NEON chunk.
  // v_inv_sigma and v_exp_sho are pre-hoisted variables to radically limit
  // computational depth loop repetition on the Cortex-A9 pipeline.
  float32x4_t calc_weight4(float32x4_t v_ox, float32x4_t v_oy,
                           float32x4_t v_pca, float32x4_t v_psa, size_t s_idx,
                           float32x4_t vmeas, float32x4_t v_inv_sigma,
                           float32x4_t v_exp_sho) {
    const auto &s = sensors[s_idx];

    // initial straight raycast to determine dynamic divergence
    // use pre-computed mounting sin/cos
    float32x4_t v_rx0 =
        vmlaq_f32(vmulq_n_f32(v_psa, s.cos_a), v_pca, vdupq_n_f32(s.sin_a));
    float32x4_t v_ry0 =
        vmlsq_f32(vmulq_n_f32(v_pca, s.cos_a), v_psa, vdupq_n_f32(s.sin_a));
    RayHit4 hit0 = raycast4(map, v_ox, v_oy, v_rx0, v_ry0);

    // evaluate the center ray (0 divergence) for the strongest-return distance
    // the peak of the Gaussian weight at the center is exactly 1.0f
    float32x4_t d0 = vminq_f32(hit0.distance, vdupq_n_f32(max_range));
    float32x4_t cos_alpha0 = vabsq_f32(vaddq_f32(
          vmulq_f32(hit0.normal_x, v_rx0), vmulq_f32(hit0.normal_y, v_ry0)));
    
    // set the best variables with the center ray before checking the rest
    float32x4_t v_max_power = vmulq_f32(cos_alpha0, fast_recip(vmulq_f32(d0, d0)));
    float32x4_t v_best_dist = d0;

    // determine dynamic mode: 24 deg if < 8 inches (mode 0), else 36 deg
    // (mode 1)
    const float32x4_t v_eight = vdupq_n_f32(8.0f);
    uint32x4_t is_near = vcltq_f32(d0, v_eight);

    // ray bundle (8 rays)
    for (int r = 0; r < 8; ++r) {
      // Pick pre-computed orientations based on proximity mode
      float32x4_t v_cm = vbslq_f32(is_near, vdupq_n_f32(s.ray_cos[0][r]),
                                   vdupq_n_f32(s.ray_cos[1][r]));
      float32x4_t v_sm = vbslq_f32(is_near, vdupq_n_f32(s.ray_sin[0][r]),
                                   vdupq_n_f32(s.ray_sin[1][r]));
      float32x4_t v_weight = vbslq_f32(is_near, vdupq_n_f32(s.ray_w[0][r]),
                                       vdupq_n_f32(s.ray_w[1][r]));

      // convert compass orientation to standard math vector
      float32x4_t v_rx = vmlaq_f32(vmulq_f32(v_psa, v_cm), v_pca, v_sm);
      float32x4_t v_ry = vmlsq_f32(vmulq_f32(v_pca, v_cm), v_psa, v_sm);

      RayHit4 hit = raycast4(map, v_ox, v_oy, v_rx, v_ry);
      float32x4_t d = vminq_f32(hit.distance, vdupq_n_f32(max_range));

      // lambertian reflection: cos_alpha = |dot(normal, -ray_dir)|
      // RayHit4 normal is already pointing towards sensor origin.
      float32x4_t cos_alpha = vabsq_f32(vaddq_f32(
          vmulq_f32(hit.normal_x, v_rx), vmulq_f32(hit.normal_y, v_ry)));

      // power model: power = (weight * cos_alpha) / d^2
      // Standard ToF power model with incidence angle
      float32x4_t v_power = vmulq_f32(vmulq_f32(v_weight, cos_alpha),
                                      fast_recip(vmulq_f32(d, d)));

      uint32x4_t better = vcgtq_f32(v_power, v_max_power);
      v_max_power = vbslq_f32(better, v_power, v_max_power);
      v_best_dist = vbslq_f32(better, d, v_best_dist);
    }

    // standard deviation stuff
    float32x4_t v_dist = v_best_dist;
    float32x4_t v_diff = vmulq_f32(vsubq_f32(vmeas, v_dist), v_inv_sigma);

    // PI_SQRT3  * (distance - max_distance) / stddev
    float32x4_t v_exp_max = fast_exp_neon(
        vmulq_f32(vsubq_f32(v_dist, vdupq_n_f32(max_range_threshold)),
                  vmulq_f32(vdupq_n_f32(PI_SQRT3), v_inv_sigma)));

    // -phi * distance
    float32x4_t v_exp_hit = fast_exp_neon(vmulq_f32(vdupq_n_f32(-phi), v_dist)),
                // -0.5 * (measurement - distance)^2 / stddev^2
        v_exp_gauss = fast_exp_neon(
            vmulq_f32(vdupq_n_f32(-0.5f), vmulq_f32(v_diff, v_diff)));

    // alpha_base = 1 - random_particle_probability
    // alpha_hit = exp_hit * alpha_base
    float32x4_t v_alpha_base = vdupq_n_f32(1.0f - random_particle_probability),
                v_alpha_hit = vmulq_f32(v_exp_hit, v_alpha_base),
                v_weight = vdupq_n_f32(0.0f);

    // distance <= max_range_threshold
    uint32x4_t m_range = vcleq_f32(v_dist, vdupq_n_f32(max_range_threshold));
    // w_hit = uniform_noise_probability * random_particle_probability *
    // alpha_hit * (1 / (sqrt(2pi) * stddev)) * exp_gauss
    float32x4_t w_hit =
        vmlaq_f32(vmulq_n_f32(vdupq_n_f32(uniform_noise_probability),
                              random_particle_probability),
                  vmulq_f32(v_alpha_hit,
                            vdupq_n_f32(INV_SQRT_2PI / distance_sensor_stddev)),
                  v_exp_gauss);

    //
    v_weight = vbslq_f32(m_range, w_hit, v_weight);

    // measurement < distance
    uint32x4_t m_short = vcltq_f32(vmeas, v_dist);
    // w_short = alpha_hit * phi * exp_short
    v_weight =
        vbslq_f32(m_short,
                  vmlaq_f32(v_weight, vmulq_f32(v_alpha_hit, vdupq_n_f32(phi)),
                            v_exp_sho),
                  v_weight);

    // measurement >= max_range_threshold
    uint32x4_t m_max = vcgeq_f32(vmeas, vdupq_n_f32(max_range_threshold));
    // w_max = alpha_hit * 1 / (1 + exp_max)
    v_weight = vbslq_f32(
        m_max,
        vaddq_f32(v_weight,
                  vmulq_f32(v_alpha_hit, fast_recip(vaddq_f32(vdupq_n_f32(1.0f),
                                                              v_exp_max)))),
        v_weight);

    return v_weight;
  }

  float wall_distance(float x0, float y0, float cos_theta, float sin_theta) {
    // convert compass direction to standard math for raycaster
    // rx = sin_theta, ry = cos_theta
    RayHit4 hit = raycast4(map, vdupq_n_f32(x0), vdupq_n_f32(y0),
                           vdupq_n_f32(sin_theta), vdupq_n_f32(cos_theta));
    return std::min(vgetq_lane_f32(hit.distance, 0), max_range);
  }

  // used to refresh the sin and cos of the angles and prevent it from being
  // larger than 1
  void refresh_trig() {
    for (size_t i = 0; i < num_aligned; i += 4) {
      // load 4 angles
      float32x4_t vt = vld1q_f32(&pt[i]);
      float32x4_t vs, vc;
      // compute sin and cos
      fast_sincos_neon(vt, vs, vc);
      // store sin and cos
      vst1q_f32(&pca[i], vc);
      vst1q_f32(&psa[i], vs);
    }
  }
};

#endif // POOPERS_ODOM_HPP_FINAL