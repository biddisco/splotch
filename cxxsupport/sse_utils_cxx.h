/*
 *  This file is part of libcxxsupport.
 *
 *  libcxxsupport is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libcxxsupport is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libcxxsupport; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libcxxsupport is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*! \file sse_utils_cxx.h
 *  SSE/SSE2/SSE3-related functionality for C++
 *
 *  Copyright (C) 2011 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_SSE_UTILS_CXX_H
#define PLANCK_SSE_UTILS_CXX_H

#include "sse_utils.h"
#include "lsconstants.h"

template<typename T, int sz> class svec;

#ifdef PLANCK_HAVE_SSE

template<> class svec<float, 4>
  {
  public:
    typedef __m128 Tv;
    typedef union { Tv v; float d[4]; } Tu;
    Tv v;

    svec () {}
    svec (const svec &b) : v(b.v) {}
    svec (const Tv &b) : v(b) {}
    svec (const float &val) : v(_mm_set1_ps(val)) {}
    svec (float val1, float val2, float val3, float val4)
      : v(_mm_set_ps(val4,val3,val2,val1)) {}
    const svec &operator= (const float &val)
      { v=_mm_set1_ps(val); return *this; }
    const svec &operator= (const svec &b)
      { v=b.v; return *this; }

    float operator[] (int p) const
      { Tu u; u.v=v; return u.d[p]; }
    void set (int p, float val)
      { Tu u; u.v=v; u.d[p]=val; v=u.v; }

    const svec &operator+= (const svec &b)
      { v=_mm_add_ps(v,b.v); return *this; }
    const svec &operator-= (const svec &b)
      { v=_mm_sub_ps(v,b.v); return *this; }
    const svec &operator*= (const svec &b)
      { v=_mm_mul_ps(v,b.v); return *this; }
    const svec &operator/= (const svec &b)
      { v=_mm_div_ps(v,b.v); return *this; }

    svec operator+ (const svec &b) const
      { return svec(_mm_add_ps(v,b.v)); }
    svec operator- (const svec &b) const
      { return svec(_mm_sub_ps(v,b.v)); }
    svec operator* (const svec &b) const
      { return svec(_mm_mul_ps(v,b.v)); }
    svec operator/ (const svec &b) const
      { return svec(_mm_div_ps(v,b.v)); }

    const svec &operator&= (const svec &b)
      { v=_mm_and_ps(v,b.v); return *this; }
    const svec &operator|= (const svec &b)
      { v=_mm_or_ps(v,b.v); return *this; }
    const svec &operator^= (const svec &b)
      { v=_mm_xor_ps(v,b.v); return *this; }
    svec operator& (const svec &b) const
      { return svec(_mm_and_ps(v,b.v)); }
    svec operator| (const svec &b) const
      { return svec(_mm_or_ps(v,b.v)); }
    svec operator^ (const svec &b) const
      { return svec(_mm_xor_ps(v,b.v)); }

    svec operator- () const
      { return svec(_mm_xor_ps(_mm_set1_ps(-0.),v)); }

    svec eq (const svec &b) const
      { return svec(_mm_cmpeq_ps(v,b.v)); }
    svec neq (const svec &b) const
      { return svec(_mm_cmpneq_ps(v,b.v)); }
    svec lt (const svec &b) const
      { return svec(_mm_cmplt_ps(v,b.v)); }
    svec le (const svec &b) const
      { return svec(_mm_cmple_ps(v,b.v)); }
    svec gt (const svec &b) const
      { return svec(_mm_cmpgt_ps(v,b.v)); }
    svec ge (const svec &b) const
      { return svec(_mm_cmpge_ps(v,b.v)); }

    void writeTo (float *val) const
      { _mm_storeu_ps (val, v); }
    void writeTo (float &a, float &b, float &c, float &d) const
      { Tu u; u.v=v; a=u.d[0]; b=u.d[1]; c=u.d[2]; d=u.d[3]; }
    void readFrom (const float *val)
      { v=_mm_loadu_ps(val); }
    void readFrom (float a, float b, float c, float d)
      { v=_mm_set_ps(d,c,b,a); }
  };

typedef svec<float,4> V4sf;

inline V4sf sqrt(const V4sf &v)
  { return V4sf(_mm_sqrt_ps(v.v)); }
inline V4sf abs(const V4sf &v)
  { return V4sf(_mm_andnot_ps(_mm_set1_ps(-0.),v.v)); }
inline V4sf blend(const V4sf &mask, const V4sf &a, const V4sf &b)
  { return V4sf(_mm_or_ps(_mm_and_ps(a.v,mask.v),_mm_andnot_ps(mask.v,b.v))); }
inline bool any (const V4sf &a)
  { return _mm_movemask_ps(a.v)!=0; }
inline bool all (const V4sf &a)
  { return _mm_movemask_ps(a.v)==15; }
inline bool none (const V4sf &a)
  { return _mm_movemask_ps(a.v)==0; }

#endif

#endif
