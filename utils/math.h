/*
 *  Ray++ - Object-oriented ray tracing library
 *  Copyright (C) 1998-2001 Martin Reinecke and others.
 *  See the AUTHORS file for more information.
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Library General Public
 *  License as published by the Free Software Foundation; either
 *  version 2 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Library General Public License for more details.
 *
 *  You should have received a copy of the GNU Library General Public
 *  License along with this library; if not, write to the Free
 *  Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  See the README file for more information.
 */

#ifndef RAYPP_MATH_H
#define RAYPP_MATH_H

#include "kernel/kernel.h"
#include "utils/twister.h"

namespace RAYPP {

extern TWISTER Rng;

inline float8 clamp (float8 val, float8 Lo, float8 Hi)
  {
  if (val < Lo) return Lo;
  if (val > Hi) return Hi;
  return val;
//  alternatively
//  return ((val<Lo) ? Lo : ((val>Hi) ? Hi : val));
//  alternatively
//  return max (Lo, min (val, Hi));
  }

inline float8 cycloidal (float8 x)
  {
  return sin (x*2*Pi);
  }

inline float8 trianglewave (float8 x)
  {
  float8 offset = fmod(x, float8(1));
  if (offset < 0)
    offset += 1;
  if (offset > 0.5)
    offset = 1 - offset;
  return offset + offset;
  }

inline float8 lerp (float8 t, float8 lo, float8 hi)
  {
  return lo + t*(hi-lo);
  }

inline float8 scurve (float8 val)
  {
  return (val*val*(3. - 2.*val));  
  }

inline float8 posmod(float8 x, float8 y)
  {
  float8 ret = fmod(x, y);
  if (ret < 0.0) ret += y;
  return ret / y;
  }

inline void seedrandom(int seed) 
  {
  srand(seed);
  }

// returns a pseudo-random float in [0,1]
inline float8 dblrand ()
  {
  return float8(rand())/RAND_MAX;
  }

// returns a pseudo-random float in [min,max]
inline float8 dblrand (float8 min, float8 max)
  {
  return (max-min) * dblrand() + min;
  }

} // namespace RAYPP

#endif
