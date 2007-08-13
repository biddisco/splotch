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

#ifndef RAYPP_AXISBOX_H
#define RAYPP_AXISBOX_H

#include "config/config.h"
#include "kernel/vector.h"
#include "kernel/interval.h"

namespace RAYPP {

class TRANSFORM;
class GEOM_RAY;

// here is still some optimization potential

/**
  \class AXISBOX kernel/axisbox.h kernel/axisbox.h
  An axis-aligned bounding box.
*/
class AXISBOX
  {
  private:
    VECTOR Min, Max;

    inline bool Check_Ray
      (const GEOM_RAY &Ray, float8 &mind, float8 &maxd) const;

  public:
    class BAD_NORMAL {};

      /*! Creates an empty box. */
    AXISBOX ();
      /*!  Creates a box with \a min as minimal coordinates and \a max
           as maximal coordinates. */
    AXISBOX (const VECTOR &min, const VECTOR &max)
      : Min (min), Max (max) {}
    AXISBOX (const IV8 &ix, const IV8 &iy, const IV8 &iz)
      : Min (ix.lo, iy.lo, iz.lo), Max (ix.hi, iy.hi, iz.hi) {}

      /*! Returns \e true if the box is empty, else \e false. */
    bool Empty () const;

      /*! Applies the linear transformation \a trans to the box.
          All eight corners of the box are transformed and the smallest
          AXISBOX enclosing them is constructed.<BR>
          \warning Don't apply more than one transformation to
          an AXISBOX, because this might unnecessarily increase the size of the
          box. It is better to gather all transformations into one and then
          apply the combined transformation. */
    void Transform (const TRANSFORM &trans);

      /*! Returns \e true if \a Loc is inside the box, else \e false. */
    bool Inside (const VECTOR &Loc) const;

      /*! Returns \e true if \a Ray lies completely or partially inside the
          box, else \e false. The region between
          \a (Ray.start+mind*Ray.dir) and \a (Ray.start+maxd*Ray.dir)
          is the intersection of the ray and the inside of the box. */
    bool Ray_in_Bounds (const GEOM_RAY &Ray, float8 &mind, float8 &maxd) const;
      /*! Same as above, but only \a mind is set. */
    bool Ray_in_Bounds (const GEOM_RAY &Ray, float8 &mind) const;
      /*! Same as above, but \a mind and \a maxd are returned in
          \a Ray.mindist and \a Ray.maxdist. */
    bool Clip_Ray (GEOM_RAY &Ray) const;

      /*! The box grows so that it includes the point \a vec. */
    void Include (const VECTOR &vec)
      { Min.Minimize (vec); Max.Maximize (vec); }
      /*! The box grows so that it includes \a box. */
    void Include (const AXISBOX &box)
      { Min.Minimize (box.Min); Max.Maximize (box.Max); }

      /*! The box is set to the intersection of itself and \a box. */
    void Build_Intersection (const AXISBOX &box)
      {
      Min.Maximize (box.Min); Max.Minimize (box.Max);
      if (Empty()) Reset();
      }

      /*! Sets the box to empty. */
    void Reset ();

      /*! Sets the minimal and maximal coordinates. */
    void Set_Corners (const VECTOR &min, const VECTOR &max)
      { Min = min; Max = max; }

      /*! Returns the total surface of the box. */
    float8 SurfaceArea () const;
      /*! Returns the maximum extent in a coordinate direction. */
    float8 MaxExtent () const
      { return max (Max.x-Min.x, max (Max.y-Min.y, Max.z-Min.z)); }

      /*! Returns the normal vector at Loc. */
    VECTOR Normal (const VECTOR &Loc) const;

      /*! Returns the minimal coordinates. */
    VECTOR Minimum () const { return Min; }
      /*! Returns the maximal coordinates. */
    VECTOR Maximum () const { return Max; }
      /*! Returns the center (average of minimal and maximal coordinates)
          of the box. */
    VECTOR Center () const { return (Min+Max)*0.5; }

      /*! Returns \e true if the box is extremely large, else \e false.*/
    bool Infinite () const;

// Union is possibly a bad name
      /*! Returns the tightest AXISBOX enclosing \a box1 and \a box2. */
    friend AXISBOX Union (const AXISBOX &box1, const AXISBOX &box2);
      /*! Returns the largest AXISBOX inside both \a box1 and \a box2. */
    friend AXISBOX Intersection (const AXISBOX &box1, const AXISBOX &box2);

      /*! Simple output operator. */
    friend ostream &operator<< (ostream &os, const AXISBOX &box);
  };

} // namespace RAYPP

#endif
