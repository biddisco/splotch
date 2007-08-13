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

#ifndef RAYPP_INTERVAL_H
#define RAYPP_INTERVAL_H

#include "config/config.h"
#include "kernel/constants.h"

namespace RAYPP {

class INTERVAL
  {
  public:
    float8 lo, hi;

    INTERVAL () {}
    INTERVAL (float8 Lo, float8 Hi) : lo(Lo), hi(Hi) {}
    INTERVAL (float8 val) : lo(val), hi(val) {}

    void Set (float8 Lo, float8 Hi)
      { lo=Lo; hi=Hi; }

    float8 len () const
      { return hi-lo; }
    float8 mid () const
      { return 0.5*(hi+lo); }
    void split (INTERVAL &i1, INTERVAL &i2) const
      {
      float8 m = mid();
      i1.Set (lo, m);
      i2.Set (m, hi);
      }

    INTERVAL operator- () const
      { return INTERVAL (-hi, -lo); }

    INTERVAL cos() const
      {
      float8 tmp1 = ::std::ceil(lo/Pi),
             tmp2 = hi/Pi,
             tmp3, tmp4;

      if (tmp1>tmp2)
        {
        if ((tmp3=::std::cos(hi)) >= (tmp4=::std::cos(lo)))
          return INTERVAL (tmp4, tmp3);
        else
          return INTERVAL (tmp3, tmp4);
        }
      else
        {
        if (tmp1+1.0 <= tmp2)
          return INTERVAL (-1.0, 1.0);
        else
          {
          if ((tmp2=tmp3=::std::cos(hi)) < (tmp4=::std::cos(lo)))
            { tmp3 = tmp4; tmp4 = tmp2; }

          if ((int4(tmp1) % 2) == 0)
            return INTERVAL (tmp4, 1.0);
          else
            return INTERVAL (-1.0, tmp3);
          }
        }
      }
    INTERVAL sin() const
      {
      return INTERVAL (lo-Pi/2, hi-Pi/2).cos();
      }
    INTERVAL atan() const
      {
      return INTERVAL (::std::atan(lo), ::std::atan(hi));
      }
    INTERVAL sinh() const
      {
      return INTERVAL (::std::sinh(lo), ::std::sinh(hi));
      }
    INTERVAL log() const
      {
      return INTERVAL (::std::log(lo), ::std::log(hi));
      }
    INTERVAL exp() const
      {
      return INTERVAL (::std::exp(lo), ::std::exp(hi));
      }
    INTERVAL sqrt () const
      {
      return INTERVAL (::std::sqrt(lo), ::std::sqrt(hi));
      }
    INTERVAL sqr () const
      {
      if (hi<0) return INTERVAL (hi*hi, lo*lo);
      if (lo>0) return INTERVAL (lo*lo, hi*hi);
      return INTERVAL (0, max(lo*lo, hi*hi));
      }

    INTERVAL &operator+= (const INTERVAL &i1)
      {
      lo += i1.lo;
      hi += i1.hi;
      return *this;
      }
    INTERVAL &operator+= (float8 f8)
      {
      lo += f8;
      hi += f8;
      return *this;
      }
    INTERVAL &operator-= (const INTERVAL &i1)
      {
      lo -= i1.hi;
      hi -= i1.lo;
      return *this;
      }
    INTERVAL &operator-= (float8 f8)
      {
      lo -= f8;
      hi -= f8;
      return *this;
      }
    INTERVAL &operator*= (const INTERVAL &i1)
      {
      float8 ac = i1.lo*lo,
             ad = i1.lo*hi,
             bc = i1.hi*lo,
             bd = i1.hi*hi;

      if (ac>ad) swap (ac, ad);
      if (bc>bd) swap (bc, bd);

      lo = min (ac, bc);
      hi = max (ad, bd);

      return *this;
      }
    INTERVAL &operator*= (float8 f8)
      {
      if (f8<0) swap (lo, hi);
      lo *= f8;
      hi *= f8;

      return *this;
      }
    INTERVAL &operator/= (const INTERVAL &i1)
      {
      if (i1.lo*i1.hi <= 0.0)
        {
        lo = -Huge_float8;
        hi =  Huge_float8;
        }
      else
        {
        (*this) *= INTERVAL (1.0/i1.hi, 1.0/i1.lo);
        }
      return *this;
      }
    INTERVAL &operator/= (float8 f8)
      {
      if (f8<0) swap (lo, hi);
      lo /= f8;
      hi /= f8;
      
      return *this;
      }
    INTERVAL operator+ (const INTERVAL &i1) const
      { return INTERVAL (lo+i1.lo, hi+i1.hi); }
    INTERVAL operator+ (float8 f8) const
      { return INTERVAL (lo+f8, hi+f8); }
    friend INTERVAL operator+ (float8 f8, const INTERVAL &i1)
      { return INTERVAL (i1.lo+f8, i1.hi+f8); }
    INTERVAL operator- (const INTERVAL &i1) const
      { return INTERVAL (lo-i1.hi, hi-i1.lo); }
    INTERVAL operator- (float8 f8) const
      { return INTERVAL (lo-f8, hi-f8); }
    friend INTERVAL operator- (float8 f8, const INTERVAL &i1)
      { return INTERVAL (f8-i1.hi, f8-i1.lo); }
    INTERVAL operator* (const INTERVAL &i1) const
      {
      INTERVAL res = *this;
      return (res*=i1);
      }
    INTERVAL operator* (float8 f8) const
      {
      INTERVAL res = *this;
      return (res*=f8);
      }
    friend INTERVAL operator* (float8 f8, const INTERVAL &i1)
      {
      INTERVAL res = i1;
      return (res*=f8);
      }
    INTERVAL operator/ (const INTERVAL &i1) const
      {
      INTERVAL res = *this;
      return (res/=i1);
      }
    INTERVAL operator/ (float8 f8) const
      {
      INTERVAL res = *this;
      return (res/=f8);
      }
    friend INTERVAL operator/ (float8 f8, const INTERVAL &i1)
      {
      if (i1.lo*i1.hi <= 0.0)
        return INTERVAL (-Huge_float8, Huge_float8);

      INTERVAL res (f8/i1.hi, f8/i1.lo);
      if (f8<0) swap (res.lo, res.hi);
      return res;
      }

    friend ostream &operator<< (ostream &os, const INTERVAL &i)
      {
      os << "[" << i.lo << "," << i.hi << "]" << endl;
      return os;
      }
  };

typedef INTERVAL IV8;

} // namespace RAYPP

#endif
