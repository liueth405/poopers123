#pragma once
/* 
 * Exponential PRNG generator.
 */
#include "./shared.h"
#define SIMDE_ENABLE_NATIVE_ALIASES
#include "simde/x86/sse2.h"

#define	__EXP_LAYERS__	252

#ifdef __cplusplus
extern "C" {
#endif

extern double __exp_X__[__EXP_LAYERS__+1];
extern uint8_t map_exp[256];
extern int64_t ipmf_exp[256];
extern double Y_exp[__EXP_LAYERS__+1];

#ifdef __cplusplus
}
#endif

static inline void exponential_setup(void){
	mt_init();
}

static inline double exponential(void); // Forward declaration

static inline double _exp_overhang(uint_fast8_t j) {
    double *X_j = __exp_X__ + j;
	MT_FLUSH();
    int64_t U_x = RANDOM_INT63();
    int64_t U_distance = RANDOM_INT63() - U_x;
    if (U_distance < 0) {
        U_distance = -U_distance;
        U_x -= U_distance;
    }
    static int64_t iE_max = 853965788476313639;
    double x = _FAST_PRNG_SAMPLE_X(X_j, U_x);   
    if (U_distance >= iE_max) return x;
    return _FAST_PRNG_SAMPLE_Y_EXP(j, pow(2, 63) - (U_x + U_distance)) <= exp(-x) ? x : _exp_overhang(j); 
}

static inline uint_fast8_t _exp_sample_A(void) {
    uint_fast8_t j = Rand->sl & 0xff;
    return Rand++->sl >= ipmf_exp[j] ? map_exp[j] : j;
}

static inline double exponential(void) {
    static uint_fast8_t i_max = __EXP_LAYERS__;
    MT_FLUSH();
    uint_fast8_t i = Rand->sl & 0xff;
    if (i < i_max) return __exp_X__[i]*RANDOM_INT63();
    Rand++;
    uint_fast8_t j = _exp_sample_A();
    static double X_0 = 7.56927469415;
    return j > 0 ? _exp_overhang(j) : X_0 + exponential();
}