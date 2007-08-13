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

#ifndef RAYPP_GEOM_RAY_H
#define RAYPP_GEOM_RAY_H

#include "config/config.h"
#include "kernel/vector.h"
#include "kernel/constants.h"

namespace RAYPP {

/**
  \class GEOM_RAY kernel/geom_ray.h kernel/geom_ray.h
  A geometrical ray in 3D space.
*/
class GEOM_RAY
  {
  public:
    /** \class DIRSTAT kernel/geom_ray.h kernel/geom_ray.h
        Helper class for GEOM_RAY. */
    class DIRSTAT
      {
      public:
        /*! */
        bool pos_x     : 1;
        /*! */
        bool pos_y     : 1;
        /*! */
        bool pos_z     : 1;
        /*! */
        bool nonzero_x : 1;
        /*! */
        bool nonzero_y : 1;
        /*! */
        bool nonzero_z : 1;
      };

    /*! */
    VECTOR start, dir;
    /*! */
    float8 mindist, maxdist;
    /*! */
    VECTOR invdir;
    /*! */
    DIRSTAT dirstat;

    /*! */
    inline void Recalc ()
      {
      if (abs (dir.x) > Small_float8)
        {
        dirstat.nonzero_x = true;
        invdir.x = 1.0/dir.x;
        dirstat.pos_x = (dir.x > 0.0);
        }
      else dirstat.nonzero_x = false;

      if (abs (dir.y) > Small_float8)
        {
        dirstat.nonzero_y = true;
        invdir.y = 1.0/dir.y;
        dirstat.pos_y = (dir.y > 0.0);
        }
      else dirstat.nonzero_y = false;

      if (abs (dir.z) > Small_float8)
        {
        dirstat.nonzero_z = true;
        invdir.z = 1.0/dir.z;
        dirstat.pos_z = (dir.z > 0.0);
        }
      else dirstat.nonzero_z = false;
      }

    /*! */
    GEOM_RAY () {}
    /*! */
    GEOM_RAY (const VECTOR &s, const VECTOR &d, float8 mind, float8 maxd)
      : start(s), dir(d), mindist(mind), maxdist(maxd)
      { Recalc(); }

    /*! */
    void Set_Direction (const VECTOR &d)
      {
      dir = d;
      Recalc();
      }

    /*! */
    VECTOR eval (const float8 dist) const
      {
      return VECTOR (start.x + dist*dir.x,
                     start.y + dist*dir.y,
                     start.z + dist*dir.z);
      }
    /*! */
    VECTOR evalmin () const
      {
      return VECTOR (start.x + mindist*dir.x,
                     start.y + mindist*dir.y,
                     start.z + mindist*dir.z);
      }
    /*! */
    VECTOR evalmax () const
      {
      return VECTOR (start.x + maxdist*dir.x,
                     start.y + maxdist*dir.y,
                     start.z + maxdist*dir.z);
      }
  };

} // namespace RAYPP

#endif
