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

#ifndef RAYPP_NOISE_H
#define RAYPP_NOISE_H

#include "kernel/kernel.h"
#include "utils/math.h"

namespace RAYPP {

class NOISE
  {
  private:
    class HASHTABLE
      {
      private:
        int2 value[256];

      public:
        HASHTABLE (uint4 seed);

        int operator() (int4 i1, int4 i2, int4 i3) const
          {
          return value[value[value[i1%256]^(i2%256)]^(i3%256)];
          }
      };

    class VECTABLE
      {
      public:
        float4 x[256], y[256], z[256], o[256];

        VECTABLE (uint4 seed);
      };

    static const HASHTABLE HTable;
    static const VECTABLE VTable;

    static const VECTOR BIGVECT;

  public:
    float8 Noise (const VECTOR &) const;
    VECTOR DNoise (const VECTOR &) const;

    float8 fBm (VECTOR, float4, float4, uint1) const;
    VECTOR DfBm (VECTOR, float4, float4, uint1) const;

    float4 Turbulence (VECTOR, float4, float4, uint1) const;
  };

class NOISE2
  {
  private:
    class HASHTABLE
      {
      private:
        uint1 value[256];

      public:
        HASHTABLE (uint4 seed);

        int operator() (int4 i1, int4 i2) const
          {
          return value[value[i1%256]^(i2%256)];
          }
      };

    class VECTABLE
      {
      public:
        float4 x[256], y[256], o[256];

        VECTABLE (uint4 seed);
      };

    static const HASHTABLE HTable;
    static const VECTABLE VTable;

  public:
    float8 Noise (float8, float8) const;

    float8 fBm (float8, float8, float4, float4, uint1) const;

    float8 Turbulence (float8, float8, float4, float4, uint1) const;
  };

} // namespace RAYPP

#endif
