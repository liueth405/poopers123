#pragma once

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
#else
#define SIMDE_ENABLE_NATIVE_ALIASES
#include "simde/arm/neon.h"
#endif

/**
 * @brief Fast exponential approximation using NEON intrinsics (Armv7-A compatible).
 *
 * @param x Input vector (4 floats).
 * @return float32x4_t Vector containing exp(x) approximations.
 *
 * @note This is a direct translation of the provided SSE logic.
 * As previously discussed, the core algorithm's step involving
 * `j = vandq_s32(t, m)` and the subsequent combination `vaddq_s32(j, ...)`
 * appears non-standard and potentially incorrect for accurately combining
 * the integer and fractional parts via IEEE 754 bit manipulation for exp.
 * Use with caution and verify its accuracy for your needs.
 */
static inline float32x4_t fast_exp_neon(float32x4_t x)
{
    float32x4_t f, p, r;
    int32x4_t t, j;

    // Constants:
    // a = (1 << 23) / log(2)
    const float32x4_t a = vdupq_n_f32(12102203.0f);
    // m = 0xff800000 (mask for sign and exponent bits)
    const int32x4_t m = vdupq_n_s32(0xff800000);
    // ttm23 = 2^(-23)
    const float32x4_t ttm23 = vdupq_n_f32(1.1920929e-7f);
    // Polynomial coefficients
    const float32x4_t c0 = vdupq_n_f32(0.3371894346f);
    const float32x4_t c1 = vdupq_n_f32(0.657636276f);
    const float32x4_t c2 = vdupq_n_f32(1.00172476f);

    // t = (int)(a * x)
    t = vcvtq_s32_f32(vmulq_f32(a, x));

    // j = t & m 
    j = vandq_s32(t, m);

    // t = t - j
    t = vsubq_s32(t, j);

    // f = ttm23 * (float)t
    f = vmulq_f32(ttm23, vcvtq_f32_s32(t));

    // Evaluate polynomial: p = (c0 * f + c1) * f + c2.
    p = vmlaq_f32(c1, c0, f);
    p = vmlaq_f32(c2, p, f);

    // Combine integer and polynomial parts via bit manipulation.
    r = vreinterpretq_f32_s32(vaddq_s32(j, vreinterpretq_s32_f32(p)));

    return r;
}
