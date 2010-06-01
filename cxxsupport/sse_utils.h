#ifndef PLANCK_SSE_UTILS_H
#define PLANCK_SSE_UTILS_H

#if (defined(__x86_64__) && defined(__SSE2__))

#define PLANCK_HAVE_SSE2

#include <xmmintrin.h>
#include <emmintrin.h>

typedef __m128 v4sf; // vector of 4 floats (SSE1)

typedef union {
  float f[4];
  v4sf  v;
} V4SF;

#endif

#endif
