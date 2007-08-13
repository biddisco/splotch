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

#ifndef RAYPP_COLOUR_H
#define RAYPP_COLOUR_H

#include "config/config.h"
#include "kernel/constants.h"

namespace RAYPP {

      /** \class COLOUR kernel/colour.h kernel/colour.h
      A class for storing RGB colour information. */
class COLOUR
  {
  public:
    /*! */
    float4 r, g, b;

    /*! */
    COLOUR () {}
    /*! */
    COLOUR (float4 rv, float4 gv, float4 bv)
      : r (rv), g (gv), b (bv) {}

    /*! */
    float4 maxcomp () const
      {
      return max (max(r,g),b);
      }

    /*! */
    float4 Intensity () const
      {
      return (0.11*b + 0.59*g + 0.30*r);
      }

    /*! */
    float4 AbsoluteSum () const
      {
      return (abs(r) + abs(g) + abs(b));
      }

    /*! */
    void Clip ()
      {
      if (g<0.0) g=0.0; else if (g>1.0) g=1.0;
      if (r<0.0) r=0.0; else if (r>1.0) r=1.0;
      if (b<0.0) b=0.0; else if (b>1.0) b=1.0;
      }

    /*! */
    bool AllSmaller (float4 thresh) const
      {
      return ((r < thresh) && (g < thresh) && (b < thresh));
      }
    /*! */
    bool AllSmaller (const COLOUR &thresh) const
      {
      return ((r < thresh.r) && (g < thresh.g) && (b < thresh.b));
      }
    /*! */
    bool TooSmall () const
      {
      return ((r < Small_float4) && (g < Small_float4) && (b < Small_float4));
      }

    /*! */
    COLOUR operator- () const
      {
      return COLOUR (-r, -g, -b);
      }

    /*! */
    COLOUR &operator*= (const COLOUR &Col)
      {
      r *= Col.r; g *= Col.g; b *= Col.b;
      return *this;
      }
    /*! */
    COLOUR &operator*= (float4 factor)
      {
      r *= factor; g *= factor; b *= factor;
      return *this;
      }
    /*! */
    COLOUR &operator+= (const COLOUR &Col)
      {
      r += Col.r; g += Col.g; b += Col.b;
      return *this;
      }
    /*! */
    COLOUR operator+ (const COLOUR &Col2) const
      {
      return COLOUR (r+Col2.r, g+Col2.g, b+Col2.b);
      }
    /*! */
    COLOUR &operator-= (const COLOUR &Col)
      {
      r -= Col.r; g -= Col.g; b -= Col.b;
      return *this;
      }
    /*! */
    COLOUR operator- (const COLOUR &Col2) const
      {
      return COLOUR (r-Col2.r, g-Col2.g, b-Col2.b);
      }
    /*! */
    COLOUR operator* (const COLOUR &Col2) const
      {
      return COLOUR (r*Col2.r, g*Col2.g, b*Col2.b);
      }
    /*! */
    COLOUR operator* (float4 factor) const
      {
      return COLOUR (r*factor, g*factor, b*factor);
      }
    /*! */
    COLOUR operator/ (float4 div) const
      {
      float4 mult = 1.0/div;
      return COLOUR (r*mult, g*mult, b*mult);
      }

    /*! */
    friend inline COLOUR operator* (float4 factor, const COLOUR &Col)
      {
      return COLOUR (Col.r*factor, Col.g*factor, Col.b*factor);
      }

    /*! */
    COLOUR exp () const
      {
      return COLOUR (::std::exp(r), ::std::exp(g), ::std::exp(b));
      }

    /*! */
    friend ostream &operator<< (ostream &os, const COLOUR &c)
      {
      os << "(" << c.r << ", " << c.g << ", " << c.b << ")";
      return os;
      }
  };

} // namespace RAYPP

#endif
