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

#ifndef RAYPP_VECTOR_H
#define RAYPP_VECTOR_H

#include "config/config.h"

namespace RAYPP {

/**
  \class VECTOR kernel/vector.h
  A 3D vector class, designed for high efficiency.
  */
class VECTOR
  {
  public:
    /*! */
    float8 x, y, z;

    /*!
      Default constructor. Components are NOT initialized!
      */
    VECTOR () {}
    /*!
      Components are set to xc, yc and zc, respectively.
      */
    VECTOR (float8 xc, float8 yc, float8 zc)
      : x(xc), y(yc), z(zc) {}

    /*!
      After this operation, the vector points in the old direction,
      but has length 1.
      \warning No check for zero-length vectors!
      */
    void Normalize ()
      {
      float8 l = 1.0/sqrt (x*x + y*y + z*z);
      x*=l; y*=l; z*=l;
      }
    /*!
      returns a vector which points in the same direction,
      but has length 1.
      \warning No check for zero-length vectors!
      */
    VECTOR Norm () const
      {
      float8 l = 1.0/sqrt (x*x + y*y + z*z);
      return VECTOR (x*l, y*l, z*l);
      }

    /*!
      returns the vector's length.
      */
    float8 Length () const
      {
      return sqrt (x*x + y*y + z*z);
      }
    /*!
      returns the square of the vector's length.
      Faster than Length()!
      */
    float8 SquaredLength () const
      {
      return (x*x + y*y + z*z);
      }

    /*! */
    VECTOR operator- () const
      {
      return VECTOR (-x, -y, -z);
      }

    /*! */
    void Flip ()
      {
      x=-x; y=-y; z=-z;
      }

    /*! */
    VECTOR &operator+= (const VECTOR &vec)
      {
      x+=vec.x; y+=vec.y; z+=vec.z;
      return *this;
      }
    /*! */
    VECTOR operator+ (const VECTOR &vec) const
      {
      return VECTOR (x+vec.x, y+vec.y, z+vec.z);
      }
    /*! */
    VECTOR &operator-= (const VECTOR &vec)
      {
      x-=vec.x; y-=vec.y; z-=vec.z;
      return *this;
      }
    /*! */
    VECTOR operator- (const VECTOR &vec) const
      {
      return VECTOR (x-vec.x, y-vec.y, z-vec.z);
      }

    /*! */
    VECTOR &operator*= (float8 factor)
      {
      x*=factor; y*=factor; z*=factor;
      return *this;
      }
    /*! */
    VECTOR operator* (float8 factor) const
      {
      return VECTOR (x*factor, y*factor, z*factor);
      }
    /*! */
    friend inline VECTOR operator* (float8 factor, const VECTOR &vec)
      {
      return VECTOR (vec.x*factor, vec.y*factor, vec.z*factor);
      }
    /*! */
    VECTOR &operator*= (const VECTOR &vec)
      {
      x*=vec.x; y*=vec.y; z*=vec.z;
      return *this;
      }
    /*! */
    VECTOR operator* (const VECTOR &vec) const
      {
      return VECTOR (x*vec.x, y*vec.y, z*vec.z);
      }
    /*! */
    VECTOR &operator/= (float8 divisor)
      {
      float8 mult = 1.0/divisor;
      x*=mult; y*=mult; z*=mult;
      return *this;
      }
    /*! */
    VECTOR operator/ (float8 divisor) const
      {
      float8 mult = 1.0/divisor;
      return VECTOR (x*mult, y*mult, z*mult);
      }
    /*! */
    VECTOR &operator/= (const VECTOR &vec)
      {
      x/=vec.x; y/=vec.y; z/=vec.z;
      return *this;
      }

    /*! */
    VECTOR operator/ (const VECTOR &vec) const
      {
      return VECTOR (x/vec.x, y/vec.y, z/vec.z);
      }
    /*! */
    void Minimize (const VECTOR &vec)
      {
      if (x > vec.x) x = vec.x;
      if (y > vec.y) y = vec.y;
      if (z > vec.z) z = vec.z;
      }
    /*! */
    void Maximize (const VECTOR &vec)
      {
      if (x < vec.x) x = vec.x;
      if (y < vec.y) y = vec.y;
      if (z < vec.z) z = vec.z;
      }

    /*! */
    friend inline float8 Dot (const VECTOR &a, const VECTOR &b)
      {
      return a.x*b.x + a.y*b.y + a.z*b.z;
      }
    /*! */
    friend inline VECTOR Cross (const VECTOR &a, const VECTOR &b)
      {
      return VECTOR (a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
      }

    /*! */
    friend ostream &operator<< (ostream &os, const VECTOR &vec)
      {
      os << "<" << vec.x << ", " << vec.y << ", " << vec.z << ">";
      return os;
      }
  };

} // namespace RAYPP

#endif
