#ifndef PLANCK_SSE_UTILS_H
#define PLANCK_SSE_UTILS_H

#if (defined(__SSE__))

#define PLANCK_HAVE_SSE

#include <xmmintrin.h>

typedef __m128 v4sf; /* vector of 4 floats (SSE1) */

typedef union {
  float f[4];
  v4sf v;
} V4SF;

static inline v4sf build_v4sf (float a, float b, float c, float d)
  {
  V4SF tmp;
  tmp.f[0]=a; tmp.f[1]=b; tmp.f[2]=c; tmp.f[3]=d;
  return tmp.v;
  }
static inline void read_v4sf (v4sf v, float *a, float *b, float *c, float *d)
  {
  V4SF tmp;
  tmp.v = v;
  if (a) *a=tmp.f[0];
  if (b) *b=tmp.f[1];
  if (c) *c=tmp.f[2];
  if (d) *d=tmp.f[3];
  }

#endif

#if (defined(__SSE2__))

#define PLANCK_HAVE_SSE2

#include <emmintrin.h>

typedef __m128d v2df; /* vector of 2 doubles (SSE2) */

typedef union {
  double d[2];
  v2df v;
} V2DF;

static inline v2df build_v2df (double a, double b)
  {
  V2DF tmp;
  tmp.d[0]=a; tmp.d[1]=b;
  return tmp.v;
  }
static inline void read_v2df (v2df v, double *a, double *b)
  {
  V2DF tmp;
  tmp.v = v;
  if (a) *a=tmp.d[0];
  if (b) *b=tmp.d[1];
  }

#endif

#endif
