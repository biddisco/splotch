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

#include "kernel/error.h"
#include "kernel/axisbox.h"
#include "kernel/transform.h"
#include "kernel/geom_ray.h"

namespace RAYPP {

inline bool AXISBOX::Check_Ray
  (const GEOM_RAY &Ray, float8 &mind, float8 &maxd) const
  {
  float8 dist;
  mind = Ray.mindist;
  maxd = Ray.maxdist;

  if (Ray.dirstat.nonzero_x)
    {
    if (Ray.dirstat.pos_x)
      {
      if ((dist = (Max.x - Ray.start.x)*Ray.invdir.x) < mind) return false;
      if (dist < maxd) maxd = dist;
      if ((dist = (Min.x - Ray.start.x)*Ray.invdir.x) > maxd) return false;
      if (dist > mind) mind = dist;
      }
    else
      {
      if ((dist = (Min.x - Ray.start.x)*Ray.invdir.x) < mind) return false;
      if (dist < maxd) maxd = dist;
      if ((dist = (Max.x - Ray.start.x)*Ray.invdir.x) > maxd) return false;
      if (dist > mind) mind = dist;
      }
    }
  else
    if ((Ray.start.x < Min.x) || (Ray.start.x > Max.x)) return false;

  if (Ray.dirstat.nonzero_y)
    {
    if (Ray.dirstat.pos_y)
      {
      if ((dist = (Max.y - Ray.start.y)*Ray.invdir.y) < mind) return false;
      if (dist < maxd) maxd = dist;
      if ((dist = (Min.y - Ray.start.y)*Ray.invdir.y) > maxd) return false;
      if (dist > mind) mind = dist;
      }
    else
      {
      if ((dist = (Min.y - Ray.start.y)*Ray.invdir.y) < mind) return false;
      if (dist < maxd) maxd = dist;
      if ((dist = (Max.y - Ray.start.y)*Ray.invdir.y) > maxd) return false;
      if (dist > mind) mind = dist;
      }
    }
  else
    if ((Ray.start.y < Min.y) || (Ray.start.y > Max.y)) return false;

  if (Ray.dirstat.nonzero_z)
    {
    if (Ray.dirstat.pos_z)
      {
      if ((dist = (Max.z - Ray.start.z)*Ray.invdir.z) < mind) return false;
      if (dist < maxd) maxd = dist;
      if ((dist = (Min.z - Ray.start.z)*Ray.invdir.z) > maxd) return false;
      if (dist > mind) mind = dist;
      }
    else
      {
      if ((dist = (Min.z - Ray.start.z)*Ray.invdir.z) < mind) return false;
      if (dist < maxd) maxd = dist;
      if ((dist = (Max.z - Ray.start.z)*Ray.invdir.z) > maxd) return false;
      if (dist > mind) mind = dist;
      }
    }
  else
    if ((Ray.start.z < Min.z) || (Ray.start.z > Max.z)) return false;

  return true;
  }


AXISBOX::AXISBOX () 
  : Min (Huge_float8, Huge_float8, Huge_float8),
    Max (-Huge_float8, -Huge_float8, -Huge_float8) {}

bool AXISBOX::Empty () const
  {
  return ((Max.x < Min.x) || (Max.y < Min.y) || (Max.z < Min.z));
  }

void AXISBOX::Transform (const TRANSFORM &trans)
  {
  VECTOR loc  = trans.TransPoint (Min),
         dirx = trans.TransDirection (VECTOR (Max.x-Min.x,0,0)),
         diry = trans.TransDirection (VECTOR (0,Max.y-Min.y,0)),
         dirz = trans.TransDirection (VECTOR (0,0,Max.z-Min.z));

  Min = Max = loc;

  loc+=dirx; Include (loc);
  loc+=diry; Include (loc);
  loc-=dirx; Include (loc);
  loc+=dirz; Include (loc);
  loc-=diry; Include (loc);
  loc+=dirx; Include (loc);
  loc+=diry; Include (loc);
  }

bool AXISBOX::Inside (const VECTOR &Loc) const
  {
  if ((Loc.x >= Min.x) && (Loc.x <= Max.x) &&
      (Loc.y >= Min.y) && (Loc.y <= Max.y) &&
      (Loc.z >= Min.z) && (Loc.z <= Max.z))
    return true;
  return false;
  } 

bool AXISBOX::Ray_in_Bounds (const GEOM_RAY &Ray, float8 &mind, float8 &maxd) const
  {
  return Check_Ray (Ray, mind, maxd);
  }

bool AXISBOX::Ray_in_Bounds (const GEOM_RAY &Ray, float8 &mind) const
  {
  float8 maxd;
  return Check_Ray (Ray, mind, maxd);
  }

bool AXISBOX::Clip_Ray (GEOM_RAY &Ray) const
  {
  float8 mind, maxd;
  if (!Check_Ray (Ray, mind, maxd)) return false;
  Ray.mindist = mind; Ray.maxdist = maxd;
  return true; 
  }

void AXISBOX::Reset ()
  {
  Min = VECTOR ( Huge_float8,  Huge_float8,  Huge_float8);
  Max = VECTOR (-Huge_float8, -Huge_float8, -Huge_float8);
  }

float8 AXISBOX::SurfaceArea () const
  {
  float8 area = (2.0 * ((Max.x-Min.x)*(Max.y-Min.y)
                      + (Max.x-Min.x)*(Max.z-Min.z)
                      + (Max.y-Min.y)*(Max.z-Min.z)));
  if (area > 0.0) return area;
  return 0.0;
  }

VECTOR AXISBOX::Normal (const VECTOR &Loc) const
  {
  if (abs(Loc.x-Min.x) < Small_dist) return VECTOR (-1,0,0);
  if (abs(Loc.x-Max.x) < Small_dist) return VECTOR ( 1,0,0);
  if (abs(Loc.y-Min.y) < Small_dist) return VECTOR (0,-1,0);
  if (abs(Loc.y-Max.y) < Small_dist) return VECTOR (0, 1,0);
  if (abs(Loc.z-Min.z) < Small_dist) return VECTOR (0,0,-1);
  if (abs(Loc.z-Max.z) < Small_dist) return VECTOR (0,0, 1);

  throw BAD_NORMAL();
  }

bool AXISBOX::Infinite () const
  {
  return ((Min.x < -Large_float8) ||
          (Max.x >  Large_float8) ||
          (Min.y < -Large_float8) ||
          (Max.y >  Large_float8) ||
          (Min.z < -Large_float8) ||
          (Max.z >  Large_float8));
  }

AXISBOX Union (const AXISBOX &box1, const AXISBOX &box2)
  {
  AXISBOX box (box1);
  box.Include (box2);
  return box;
  }

AXISBOX Intersection (const AXISBOX &box1, const AXISBOX &box2)
  {
  AXISBOX box (box1);
  box.Build_Intersection (box2);
  return box;
  }

ostream &operator<< (ostream &os, const AXISBOX &box)
  {
  os << "Axisbox {" << box.Min << ", " << box.Max << "}";
  return os;
  }

} // namespace RAYPP
