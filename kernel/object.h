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

#ifndef RAYPP_OBJECT_H
#define RAYPP_OBJECT_H

#include "config/config.h"
#include "kernel/initable.h"
#include "kernel/transformable.h"
#include "kernel/axisbox.h"

namespace RAYPP {

class OBJECT_QUEUE;
class RAY;
class VECTOR;
class INSIDE_INFO;
class INTERSECT_INFO;

/**
  \class OBJECT kernel/object.h
  The building block of the scene.
*/
class OBJECT: public INITABLE, public TRANSFORMABLE
  {
  public:
    /*! */
    virtual AXISBOX BBox () const = 0;

    /*! */
    virtual bool Has_Volume () const = 0;
    /*! */
    virtual bool Volume_in_BBox () const = 0;

    /*! */
    virtual bool Test (RAY &Ray, float8 &mindist) const = 0;

    /*! */
    virtual bool Inside_Volume (const VECTOR &Loc,
      const INSIDE_INFO &nowInside, INSIDE_INFO &result) const = 0;

    /*! */
    virtual bool Intersect (RAY &Ray, OBJECT_QUEUE &Queue,
      INTERSECT_INFO &result) const = 0;
  };

} // namespace RAYPP

#endif
