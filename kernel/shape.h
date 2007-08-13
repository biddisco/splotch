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

#ifndef RAYPP_SHAPE_H
#define RAYPP_SHAPE_H

#include "config/config.h"
#include "kernel/initable.h"
#include "kernel/transformable.h"
#include "kernel/axisbox.h"

namespace RAYPP {

class VECTOR;
class GEOM_RAY;

/**
  \class SHAPE kernel/shape.h kernel/shape.h
  Prototype for all geometrical shapes.
  This class contains functions for determining the intersection(s) of
  a ray with a user-defined geometrical object and functions for testing
  if a 3D point is inside this object.
 */
class SHAPE: public INITABLE, public TRANSFORMABLE
  {
  public:
    /*! */
    typedef pair<float8, VECTOR> INTER;

    /*!
      Must return an AXISBOX that encloses the object completely (the tighter
      the box is, the faster the ray tracer will be). <br>
      <strong>Note:</strong> It is safe to assume that BBox() will not
      be called very often
      during image generation, and is therefore not time-critical. The idea
      is to calculate the bounding box every time BBox() is called and not
      to store it in the SHAPE. This can save significant amounts of memory.

      May only be called after initialization.
     */
    virtual AXISBOX BBox () const = 0;

    /*!
      Must return \e true, if the object divides 3D space into an inside and
      an outside region (like a sphere). <br>
      Must return \e false if the object doesn't (like a triangle).

      May only be called after initialization.
     */
    virtual bool Has_Inside () const = 0;
    /*!
      Must return \e true, if the whole inside region of the object is also
      inside the object's bounding box, or if the object doesn't have an inside
      region (regular sphere, triangle). <br>
      Must return \e false if the bounding box doesn't enclose all of the
      inside region. For example, this would be the case for an inverted sphere
      (a sphere turned 'inside-out').

      May only be called after initialization.
     */
    virtual bool Inside_in_BBox () const = 0;

    /*!
      This one is a bit complicated, but does a great job to speed up the
      ray tracing.

      Perform a quick test if Ray can hit the object at all.

      <ul><li> If not, return \e false. </li>
      <li>If yes, set mindist to a lower limit of the depth of the first
          intersection.
          If \b mindist is the exact intersection depth, set \b realhit
          to \e true, else to \e false. </li>
      </ul>

      Remember: \b mindist must lie between \b Ray.mindist
      and  <strong>Ray.maxdist</strong>!

      May only be called after initialization.
     */
    virtual bool Test (const GEOM_RAY &Ray, float8 &mindist, bool &realhit)
      const = 0;
// oder aber float8 &realdist; >Ray.maxdist, wenn kein realhit;

    /*!
      If the object has no inside, return \e false. <br>
      Return \e true, if \b Loc is inside the object, else \e false.

      May only be called after initialization.
     */
    virtual bool Inside (const VECTOR &Loc) const = 0;

// evtl. u,v
    /*!
      If \b Ray intersects the object (between \b Ray.mindist
      and <strong>Ray.maxdist</strong>), set \b dist to the depth of the first
      intersection and \b Normal to the normal direction of the surface
      at that point (\b Normal must be of unit length) and return \e true.

      If not, return false.

      If the shape has an inside, the normal should point to the outside;
      else it should point consistently to one side.

      May only be called after initialization.
     */
    virtual bool Intersect (const GEOM_RAY &Ray, float8 &dist,
      VECTOR &Normal) const = 0;

// evtl. u,v
    /*!
      Insert the depths of all intersections (between \b Ray.mindist
      and <strong>Ray.maxdist</strong>) and the normal directions
      of the surface at the intersection points into <strong>Inter</strong>.
      The normals must be of unit length.
      The normal should point consistently to one side.
      If the shape has an inside, the normals should point to the outside;
      else they should point consistently to one side.
     */
    virtual void All_Intersections (const GEOM_RAY & Ray,
      vector<INTER> &Inter) const = 0;
  };

} // namespace RAYPP

#endif
