/*
 *  Ray++ - Object-oriented ray tracing library
 *  Copyright (C) 1998-2004 Martin Reinecke and others.
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

#ifndef RAYPP_TWISTER_H
#define RAYPP_TWISTER_H

#include "config/config.h"

namespace RAYPP {

class TWISTER
  {
  private:
    typedef vector<uint4> State;
    typedef State::iterator Iter;

    class BitMixer
      {
      public:
        enum { K = 0x9908b0df };
        uint4 s0;

        BitMixer()
          : s0(0) {}

        inline uint4 low_mask (uint4 s1) 
          { return (s1&1u) ? K : 0u; }
        inline uint4 high_mask (uint4 s1) const
          { return ((s0&0x80000000)|(s1&0x7fffffff)) >> 1; }
        inline uint4 operator() (uint4 s1)
          {
          uint4 r = high_mask(s1) ^ low_mask(s1);
          s0 = s1;
          return r;
          }
      };

    enum { N = 624, PF = 397, reference_seed = 4357 };
    State       S;
    Iter        I;

  public:
    TWISTER() {}

    void seed (uint4 seed_ = reference_seed)
      {
      const uint4 Knuth_A = 69069;
      if (!S.size()) S.resize(N);
      uint4 x = seed_ & 0xFFFFFFFF;
      Iter s = S.begin();
      uint4 mask = (seed_ == reference_seed) ? 0 : 0xFFFFFFFF;
      for (uint4 j = 0; j < N; ++j)
        {
        *s++ = (x + (mask & j)) & 0xFFFFFFFF; 
        x *= Knuth_A;
        }

      I = S.begin();
      }

    void reload ()
      {
      if (!S.size()) seed (); // auto-seed detection

      Iter p0 = S.begin();
      Iter pM = p0 + PF;
      BitMixer twist;
      twist (S[0]); // prime the pump
      for (Iter pf_end = S.begin()+N-PF; p0 != pf_end; ++p0, ++pM)
        *p0 = *pM ^ twist (p0[1]);
      pM = S.begin();
      for (Iter s_end = S.begin()+N-1; p0 != s_end; ++p0, ++pM)
        *p0 = *pM ^ twist (p0[1]);
      *p0 = *pM ^ twist (S[0]);

      I = S.begin();
      }

    uint4 u4rand ()
      {
      if (I >= S.end()) reload();
      uint4 y = *I++;
      y ^= (y >> 11);
      y ^= (y <<  7) & 0x9D2C5680;
      y ^= (y << 15) & 0xEFC60000;
      y ^= (y >> 18);
      return y;
      }

    uint4 u4rand (uint4 range)
      {
      return u4rand() % range;
      }

    uint4 operator()(uint4 range)
      {
      return u4rand() % range;
      }

    float8 f8rand ()
      {
      const float8 f = (1.0 / 65536) / 65536;
      return f*(u4rand() + f*u4rand());
      }

    float8 f8rand (float8 lo, float8 hi)
      {
      return lo + (hi-lo)*f8rand();
      }
  };

} // namespace RAYPP

#endif
