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

#ifndef RAYPP_TRANSFORMABLE_H
#define RAYPP_TRANSFORMABLE_H

#include "config/config.h"
#include "kernel/transform.h"

namespace RAYPP {

class VECTOR;

/**
  \class TRANSFORMABLE kernel/transformable.h kernel/transformable.h
  This class is inherited by all classes that describe things in 3D-space.
  It allows to apply any linear 3D transformation to an object. Degenerate
  transformations must not be used.
*/
class TRANSFORMABLE
  {
  public:
    /*!
       virtual destructor to ensure proper deallocation of dynamically
       allocated objects
      */
    virtual ~TRANSFORMABLE () {}

    /*!
       This function must apply the linear transformation trans to
       the object.
      */
    virtual void Transform (const TRANSFORM &trans) = 0;

    /*!
       Frontend for Transform().
       Translates the object by <em>vec</em>.
      */
    void Translate (const VECTOR &vec)
      { Transform (Translation_Transform (vec)); }
    /*!
       Frontend for Transform().
       Rotates the object first <em>vec.x</em>
       degrees around the x-axis, then <em>vec.y</em> degrees around the
       y-axis, then <em>vec.z</em> degrees around the z-axis
       (left-handed rotation).
      */
    void Rotate (const VECTOR &vec)
      { Transform (Rotation_Transform (vec)); }
    /*!
       Frontend for Transform().
       Scales the object by <em>vec</em>.
       The vector elements must be positive.
      */
    void Scale (const VECTOR &vec)
      { Transform (Scaling_Transform (vec)); }
    /*!
       Frontend for Transform().
       Scales the object by <em>val</em>.
       <em>val</em> must be positive.
      */
    void Scale (float8 val)
      { Transform (Scaling_Transform (VECTOR (val, val, val))); }
    /*!
       Frontend for Transform().
       Rotates the object <em>angle</em>
       degrees around the vector <em>axis</em>.
      */
    void Axis_Rotate (const VECTOR &axis, float8 angle)
      { Transform (Axis_Rotation_Transform (axis, angle)); }
    /*!
       Frontend for Transform().
       Performs a shearing transformation on the object.
       The transformation looks as follows:

       \code
       xnew = x + xy*y + xz*z
       ynew = yx*x + y + yz*z
       znew = zx*x + zy*y + z 
       \endcode
      */
    void Shear
      (float4 xy, float4 xz, float4 yx, float4 yz, float4 zx, float4 zy)
      { Transform (Shearing_Transform (xy, xz, yx, yz, zx, zy)); }
  };

} // namespace RAYPP

#endif
