#pragma once
/* 
 * Normal PRN generator. Must call normal_setup() to initialize the PRNG.
 */
#include "./shared.h"
#include "exponential.h"
#define SIMDE_ENABLE_NATIVE_ALIASES
#include "simde/x86/sse2.h"
#define __NORM_BINS__	253

static inline void normal_setup(void) {
  mt_init();
}

#ifdef __cplusplus
extern "C" {
#endif

extern uint8_t map_normal[256];
extern int64_t ipmf_normal[256];
extern double X_normal[__NORM_BINS__+1];
extern double Y_normal[__NORM_BINS__+1];

#ifdef __cplusplus
}
#endif

static inline uint_fast8_t _norm_sample_A(void) {
    uint_fast8_t j = Rand->sl & 0xff;
    return Rand++->sl >= ipmf_normal[j] ? map_normal[j] : j;
}

static inline double normal(void) {
    static uint_fast16_t i_max = __NORM_BINS__;
    MT_FLUSH();
	uint_fast8_t i = Rand->l & 0xff;
	if (i < i_max) return X_normal[i]*Rand++->sl;
    uint64_t U_1 = RANDOM_INT63();
    double sign_bit = U_1 & 0x100 ? 1. : -1.;
	uint_fast8_t j = _norm_sample_A();
    int64_t U_diff;
    static int32_t max_iE = 528335327, min_iE =  177059254;
	static double X_0 = 3.6360066255;
    static uint_fast8_t j_inflection = 204;
    double x, *X_j = X_normal + j;
    
    if (j > j_inflection) {
        for (;;) {
            x = _FAST_PRNG_SAMPLE_X(X_j, U_1);
            MT_FLUSH();
            U_diff = RANDOM_INT63() - U_1;
            if (U_diff >= 0) break;      
            if (U_diff >= -max_iE && _FAST_PRNG_SAMPLE_Y_NORMAL(j, pow(2, 63) - (U_1 + U_diff)) < exp(-0.5*x*x) ) break;
            U_1 = RANDOM_INT63();
        }
    } else if (j == 0) {
        do x = pow(X_0, -1)*exponential();
        while (exponential() < 0.5*x*x);
        x += X_0;
    } else if (j < j_inflection) {
        for (;;) {
            MT_FLUSH();
            U_diff = RANDOM_INT63() - U_1;
            if (U_diff < 0) { U_diff = -U_diff; U_1 -= U_diff; }
            x = _FAST_PRNG_SAMPLE_X(X_j, U_1);
            if (U_diff > min_iE) break;
            if (_FAST_PRNG_SAMPLE_Y_NORMAL(j, pow(2, 63) - (U_1 + U_diff)) < exp(-0.5*x*x) ) break;
            U_1 = RANDOM_INT63();
        } 
    } else {
        for (;;) {
            x = _FAST_PRNG_SAMPLE_X(X_j, U_1);
            MT_FLUSH();
            if (_FAST_PRNG_SAMPLE_Y_NORMAL(j, RANDOM_INT63()) < exp(-0.5*x*x) ) break;
            U_1 = RANDOM_INT63();
        }
    }
    return sign_bit*x; 
}