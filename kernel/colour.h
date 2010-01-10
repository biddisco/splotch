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

/** \class COLOUR kernel/colour.h kernel/colour.h
    A class for storing RGB colour information. */
class COLOUR
  {
  public:
    /*! */
    float32 r, g, b;

    /*! */
    COLOUR () {}
    /*! */
    COLOUR (float32 rv, float32 gv, float32 bv)
      : r (rv), g (gv), b (bv) {}

    /*! */
    COLOUR operator+ (const COLOUR &Col2) const
      { return COLOUR (r+Col2.r, g+Col2.g, b+Col2.b); }
    /*! */
    COLOUR operator- (const COLOUR &Col2) const
      { return COLOUR (r-Col2.r, g-Col2.g, b-Col2.b); }
    /*! */
    COLOUR operator* (float32 factor) const
      { return COLOUR (r*factor, g*factor, b*factor); }
    /*! */
    friend inline COLOUR operator* (float32 factor, const COLOUR &Col)
      { return COLOUR (Col.r*factor, Col.g*factor, Col.b*factor); }

    /*! */
    friend std::ostream &operator<< (std::ostream &os, const COLOUR &c)
      {
      os << "(" << c.r << ", " << c.g << ", " << c.b << ")";
      return os;
      }
  };

#endif
