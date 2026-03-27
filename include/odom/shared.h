#pragma once
#include <stdlib.h>
#include <math.h>
#include "MT19937.h"
#define SIMDE_ENABLE_NATIVE_ALIASES
#include "simde/x86/sse2.h"

#define MAX_INT64   0x7fffffffffffffff
#define RANDOM_INT63() ( Rand++->sl & MAX_INT64 )

#define _FAST_PRNG_SAMPLE_X(X_j, U) ( *(X_j)*pow(2, 63) + ((X_j)[-1] - *(X_j) )*(U) )

#define _FAST_PRNG_SAMPLE_Y_NORMAL(i, U) ( Y_normal[(i)-1]*pow(2, 63) + (Y_normal[(i)  ] - Y_normal[(i)-1])*(U) )
#define _FAST_PRNG_SAMPLE_Y_EXP(i, U) ( Y_exp[(i)-1]*pow(2, 63) + (Y_exp[(i)  ] - Y_exp[(i)-1])*(U) )
