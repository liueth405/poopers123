#pragma once
#include <unistd.h>
#include <time.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#define SIMDE_ENABLE_NATIVE_ALIASES
#include "simde/x86/sse2.h"

#define __cycle__		500	
#define MEXP	19937
#define iN (MEXP / 128 + 1)
#define N32 (iN * 4)
#define N64 (iN * 2)

#define POS1	122
#define SL1	18
#define SL2	1
#define SR1	11
#define SR2	1
#define MSK1	0xdfffffefU
#define MSK2	0xddfecb7fU
#define MSK3	0xbffaffffU
#define MSK4	0xbffffff6U
#define PARITY1	0x00000001U
#define PARITY2	0x00000000U
#define PARITY3	0x00000000U
#define PARITY4	0x13c9e684U

#ifndef PRE_ALWAYS
#define PRE_ALWAYS inline
#endif

union W128_T {
    __m128i si;
    uint32_t u[4];
    double d[2];
};
typedef union W128_T w128_t;

typedef union RAND64_t {
	uint8_t  s[8];
	uint32_t i[2];
	double d;
	float f[2];
	uint64_t l;
	signed long long sl;
} rand64_t;

typedef union dW128_T {
    __m128i si;
    __m128d sd;
    uint64_t l[2];
    uint32_t i[4];
    double d[2];
} dw128_t;

#ifdef __cplusplus
extern "C" {
#endif

/* Global state (Moved to src/cold_math.cpp) */
extern w128_t sfmt[iN];
extern uint32_t *psfmt32;
extern int idx;
extern uint32_t parity[4];

void gen_rand_array(w128_t *array, int size);
void mt_init(void);

#define __EXP_SET__		0x3ff0000000000000
extern w128_t iRandS[__cycle__];
extern w128_t *iRend; 
extern rand64_t *Rand;
extern __m128d sse2_double_m_one;
extern __m128i sse2_int_set;

#ifdef __cplusplus
}
#endif

#ifdef REPORT_PRNS
extern long __n_cycles__;
void _report_PRN_total(void);
#define INCREMENT_N_CYCLES() (__n_cycles__++);
#else
#define INCREMENT_N_CYCLES() ;
#endif

#define MT_FLUSH() { if (Rand > (rand64_t *)iRend) { \
gen_rand_array(iRandS,__cycle__); \
Rand = (rand64_t *) iRandS; \
INCREMENT_N_CYCLES() \
}; } 

static inline dw128_t wide_uniform(void) {
  MT_FLUSH();
  dw128_t W;
  W.si = _mm_set_epi64x(Rand[0].l, Rand[1].l);	
  Rand+=2;
  W.si = _mm_or_si128(_mm_srli_epi64(W.si, 2), sse2_int_set); 
  W.sd = _mm_add_pd(W.sd, sse2_double_m_one);
  return W;
}

static inline double uniform_double_PRN(void) {
  MT_FLUSH();
  Rand->l = (Rand->l >> 2) | __EXP_SET__;
  return Rand++->d - 1;
}

static inline unsigned long long rand_long(unsigned long long n){
  MT_FLUSH();
  return Rand++->l % n;
}

static inline unsigned long long rand_long64(void){
  MT_FLUSH();
  return Rand++->l;
}
