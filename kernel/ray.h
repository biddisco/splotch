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

#ifndef RAYPP_RAY_H
#define RAYPP_RAY_H

#include "config/config.h"
#include "kernel/geom_ray.h"
#include "kernel/vector.h"
#include "kernel/colour.h"
#include "kernel/inside_info.h"

namespace RAYPP {

/**
  \class RAY kernel/ray.h
  Child of GEOM_RAY, with some additional data.
*/
class RAY: public GEOM_RAY
  {
  public:
// Inside is valid at the location start+mindist*dir
    /*! */
    INSIDE_INFO Inside;    
    /*! */
    COLOUR Importance;
    /*! */
    uint1 diffuse_level, specular_level;

    /*! */
    RAY () {}
    /*! */
    RAY (const VECTOR &s, const VECTOR &d, float8 mind, float8 maxd,
         const INSIDE_INFO &ii, const COLOUR &c, uint1 diff=0, uint1 spec=0)
      :  GEOM_RAY (s, d, mind, maxd), Inside (ii), Importance (c),
         diffuse_level (diff), specular_level (spec)
      {}

    /*! */
    void SetAll (const VECTOR &s, const VECTOR &d, float8 mind, float8 maxd,
         const INSIDE_INFO &ii, const COLOUR &c, uint1 diff=0, uint1 spec=0)
      {
      start = s;
      Set_Direction (d);
      mindist = mind;
      maxdist = maxd;
      Inside = ii;
      Importance = c;
      diffuse_level = diff;
      specular_level = spec;
      }
  };

} // namespace RAYPP

#endif
