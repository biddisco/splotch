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

#include "utils/noise.h"

namespace RAYPP {

NOISE::HASHTABLE::HASHTABLE (uint4 seed)
  {
  for (int2 i=0; i<256; ++i) value[i] = i;
  Rng.seed (seed);
  random_shuffle (&value[0], &value[0]+256, Rng);
  }

const NOISE::HASHTABLE NOISE::HTable (14);

NOISE::VECTABLE::VECTABLE (uint4 seed)
  {
  float8 l;
  Rng.seed (seed);
  for (int2 m=0; m<256; ++m)
    {
    do
      {
      x[m] = Rng.f8rand (-1, 1);
      y[m] = Rng.f8rand (-1, 1);
      z[m] = Rng.f8rand (-1, 1);
      o[m] = Rng.f8rand (-0.5, 0.5);
      }
    while ((l = x[m]*x[m] + y[m]*y[m] + z[m]*z[m]) > 1.0);
    l = 1.0/sqrt (l);
    x[m] *= l; y[m] *= l; z[m] *= l;
    }
  }

const NOISE::VECTABLE NOISE::VTable (14);

const VECTOR NOISE::BIGVECT (1000.2, 1000.5, 1000.7);

float8 NOISE::Noise (const VECTOR &vec) const
  {
  float8 x=vec.x+1e6, y=vec.y+1e6, z=vec.z+1e6;
  int ix=int(x), iy=int(y), iz=int(z);

  float8 fx0=x-ix, fx1=fx0-1, fy0=y-iy, fy1=fy0-1, fz0=z-iz, fz1=fz0-1;
  float8 wx = scurve(fx0), wy = scurve(fy0), wz = scurve(fz0);

  float8 vx0, vx1, vy0, vy1, vz0, vz1;
  int i;

  i = HTable (ix, iy, iz);
  vx0 = VTable.x[i]*fx0 + VTable.y[i]*fy0 + VTable.z[i]*fz0 + VTable.o[i];
  i = HTable (ix+1, iy, iz);
  vx1 = VTable.x[i]*fx1 + VTable.y[i]*fy0 + VTable.z[i]*fz0 + VTable.o[i];
  vy0 = vx0 + wx*(vx1-vx0);

  i = HTable (ix, iy+1, iz);
  vx0 = VTable.x[i]*fx0 + VTable.y[i]*fy1 + VTable.z[i]*fz0 + VTable.o[i];
  i = HTable (ix+1, iy+1, iz);
  vx1 = VTable.x[i]*fx1 + VTable.y[i]*fy1 + VTable.z[i]*fz0 + VTable.o[i];
  vy1 = vx0 + wx*(vx1-vx0);
  vz0 = vy0 + wy*(vy1-vy0);

  i = HTable (ix, iy, iz+1);
  vx0 = VTable.x[i]*fx0 + VTable.y[i]*fy0 + VTable.z[i]*fz1 + VTable.o[i];
  i = HTable (ix+1, iy, iz+1);
  vx1 = VTable.x[i]*fx1 + VTable.y[i]*fy0 + VTable.z[i]*fz1 + VTable.o[i];
  vy0 = vx0 + wx*(vx1-vx0);

  i = HTable (ix, iy+1, iz+1);
  vx0 = VTable.x[i]*fx0 + VTable.y[i]*fy1 + VTable.z[i]*fz1 + VTable.o[i];
  i = HTable (ix+1, iy+1, iz+1);
  vx1 = VTable.x[i]*fx1 + VTable.y[i]*fy1 + VTable.z[i]*fz1 + VTable.o[i];
  vy1 = vx0 + wx*(vx1-vx0);
  vz1 = vy0 + wy*(vy1-vy0);

  return (vz0 + wz*(vz1-vz0));
  }

VECTOR NOISE::DNoise (const VECTOR &vec) const
  {
  return VECTOR (Noise (vec-BIGVECT), Noise (vec), Noise (vec+BIGVECT));
  }

float8 NOISE::fBm (VECTOR loc, float4 lambda, float4 omega,
  uint1 Octaves) const
  {
  float4 factor=1;
  float8 val = Noise (loc); 
  for (int i=0; i<(Octaves-1); ++i)
    {
    loc *= lambda;
    factor *= omega;
    val += factor*Noise (loc);
    }
  return val;
  }

VECTOR NOISE::DfBm (VECTOR loc, float4 lambda, float4 omega,
  uint1 Octaves) const
  {
  float4 factor=1;
  VECTOR val = DNoise (loc); 
  for (int i=0; i<(Octaves-1); ++i)
    {
    loc *= lambda;
    factor *= omega;
    val += factor*DNoise (loc);
    }
  return val;
  }

float4 NOISE::Turbulence (VECTOR loc, float4 lambda, float4 omega,
  uint1 Octaves) const
  {
  float4 factor=1;
  float8 val = abs(Noise (loc)); 
  for (int i=0; i<(Octaves-1); ++i)
    {
    loc *= lambda;
    factor *= omega;
    val += factor*abs(Noise (loc));
    }
  return val;
  }


NOISE2::HASHTABLE::HASHTABLE (uint4 seed)
  {
  for (int2 i=0; i<256; ++i) value[i] = i;
  Rng.seed (seed);
  random_shuffle (&value[0], &value[0]+256, Rng);
  }

const NOISE2::HASHTABLE NOISE2::HTable (14);

NOISE2::VECTABLE::VECTABLE (uint4 seed)
  {
  float8 l;
  Rng.seed (seed);
  for (int2 m=0; m<256; ++m)
    {
    do
      {
      x[m] = Rng.f8rand (-1, 1);
      y[m] = Rng.f8rand (-1, 1);
      o[m] = Rng.f8rand (-0.5, 0.5);
      }
    while ((l = x[m]*x[m] + y[m]*y[m]) > 1.0);
    l = 1.0 / sqrt (l);
    x[m] *= l; y[m] *= l;
    }
  }

const NOISE2::VECTABLE NOISE2::VTable (14);

float8 NOISE2::Noise (float8 x, float8 y) const
  {
  x+=1e6;
  y+=1e6;
  int ix=int(x), iy=int(y);

  float8 fx0=x-ix, fx1=fx0-1, fy0=y-iy, fy1=fy0-1;
  float8 wx = scurve(fx0), wy = scurve(fy0);

  float8 vx0, vx1, vy0, vy1;
  int i;

  i = HTable (ix, iy);
  vx0 = VTable.x[i]*fx0 + VTable.y[i]*fy0 + VTable.o[i];
  i = HTable (ix+1, iy);
  vx1 = VTable.x[i]*fx1 + VTable.y[i]*fy0 + VTable.o[i];
  vy0 = vx0 + wx*(vx1-vx0);

  i = HTable (ix, iy+1);
  vx0 = VTable.x[i]*fx0 + VTable.y[i]*fy1 + VTable.o[i];
  i = HTable (ix+1, iy+1);
  vx1 = VTable.x[i]*fx1 + VTable.y[i]*fy1 + VTable.o[i];
  vy1 = vx0 + wx*(vx1-vx0);

  return (vy0 + wy*(vy1-vy0));
  }

float8 NOISE2::fBm (float8 x, float8 y, float4 lambda, float4 omega,
  uint1 Octaves) const
  {
  float8 val = Noise (x, y); 
  float4 factor = 1;
  for (int i=0; i<(Octaves-1); ++i)
    {
    x *= lambda;
    y *= lambda;
    factor *= omega;
    val += factor*Noise (x, y);
    }
  return val;
  }

float8 NOISE2::Turbulence (float8 x, float8 y, float4 lambda, float4 omega,
  uint1 Octaves) const
  {
  float4 factor = 1;
  float8 val = abs(Noise (x, y)); 
  for (int i=0; i<(Octaves-1); ++i)
    {
    x *= lambda;
    y *= lambda;
    factor *= omega;
    val += factor*abs(Noise (x, y));
    }
  return val;
  }

} // namespace RAYPP
