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

#ifndef RAYPP_VECTOR_MATH_H
#define RAYPP_VECTOR_MATH_H

#include "kernel/kernel.h"
#include "utils/math.h"

namespace RAYPP {

inline void Calc_Orthogonal_System (VECTOR v, VECTOR &v1, VECTOR &v2)
  {
  float8 l;
  if ((l = v.y*v.y + v.z*v.z) > Small_float8)
    {
    l = 1.0/sqrt(l);
    v1.x = 0.0; v1.y = v.z*l; v1.z = -v.y*l;
    v2.x = v.y*v1.z - v.z*v1.y; v2.y = - v.x*v1.z; v2.z = v.x*v1.y;
    }
  else
    {
    l = 1.0/sqrt(v.x*v.x + v.z*v.z);
    v1.x = -v.z*l; v1.y = 0.0; v1.z = v.x*l;
    v2.x = v.y*v1.z; v2.y = v.z*v1.x - v.x*v1.z; v2.z = - v.y*v1.x;
    }
  }

/*!
   Takes vector \a v and its orthognonal vectors \a ortho1 and \a ortho2 and
   returns a vector you get by rotating \a v by \a theta and \a phi.
*/
inline VECTOR Change_Dir (VECTOR v, VECTOR ortho1, VECTOR ortho2, 
  float8 theta, float8 phi) 
  {
  float8 x,y,z;
  VECTOR a (cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta));

  x = ortho1.x * a.x + ortho2.x * a.y + v.x * a.z;
  y = ortho1.y * a.x + ortho2.y * a.y + v.y * a.z;
  z = ortho1.z * a.x + ortho2.z * a.y + v.z * a.z;

  return VECTOR(x,y,z);
  }

/*!
   Calculates the direction of the reflected ray, with \a Incident_Dir
   being the incoming direction and \a Norm the surface normal. <BR>
   <tt>cosine = - Dot (Normal, Incident_Dir)</tt>
*/
inline void Calc_Reflected_Dir (VECTOR Incident_Dir, VECTOR Normal, 
                                float8 cosine, VECTOR &Reflected)
  {
  Reflected = Incident_Dir + 2 * Normal * cosine;
  }

/*!
   Calculates the direction of the reflected ray, with \a Incident_Dir
   being the incoming direction and \a Norm the surface normal.
*/
inline void Calc_Reflected_Dir (VECTOR Incident_Dir, VECTOR Normal, 
                                VECTOR &Reflected)
  {
  Reflected = Incident_Dir - 2 * Normal *  Dot(Normal, Incident_Dir);
  }

/*!
   Calculates the direction of the refracted ray, with \a Incident_Dir
   being the incoming direction, \a Norm the surface normal and \a n12
   the ratio of the inner refraction index to the outer one.
   In case of a total inner reflection this function returns false
   and Refracted is left unchanged. <BR>
   <TT>cosine = - Dot (Normal, Incident_Dir)</TT>
*/
inline bool Calc_Refracted_Dir (VECTOR Incident_Dir, VECTOR Normal, 
                                float8 cosine, float8 n12, VECTOR &Refracted)
  {
  if (abs(n12 - 1.0) > Small_float8)
    {
    // calculate refracted ray direction
    // formula from Craig Lindley, Practical Ray Tracing in C
    float8 discrim = 1 + n12*n12 * (cosine*cosine - 1);

    if (discrim<0.0)  // total inner reflection
      return false;
    else
      {
      Refracted =
        (n12 * Incident_Dir + (n12 * cosine - sqrt(discrim)) * Normal).Norm();
      }
    }
  else Refracted = Incident_Dir;

  return true;
  }

/*!
   Calculates the direction of the refracted ray, with \a Incident_Dir
   being the incoming direction, \a Norm the surface normal and \a n12
   the ratio of the inner refraction index to the outer one.
   In case of a total inner reflection this function returns false
   and \a Refracted is left unchanged.
*/
inline bool Calc_Refracted_Dir (VECTOR Incident_Dir, VECTOR Normal, 
                                float8 n12, VECTOR &Refracted)
  {
  float8 cosine = -Dot (Normal, Incident_Dir);
  if (abs(n12 - 1.0) > Small_float8)
    {
    // calculate refracted ray direction
    // formula from Craig Lindley, Practical Ray Tracing in C
    float8 discrim = 1 + n12*n12 * (cosine*cosine - 1);

    if (discrim<0.0)  // total inner reflection
      return false;
    else
      {
      Refracted =
        (n12 * Incident_Dir + (n12 * cosine - sqrt(discrim)) * Normal).Norm();
      }
    }
  else Refracted = Incident_Dir;

  return true;
  }

inline bool Intersect_Quadrangle (const VECTOR &p1, const VECTOR &p2,
  const VECTOR &p3, const VECTOR &p4, const GEOM_RAY &Ray, float8 &dist,
  float8 &ru, float8 &rv)
  {
  bool found = false;

  dist = Ray.maxdist;

  VECTOR d1 = p2 - p1,
         d2 = p3 - p1;
  VECTOR PVec = Cross (Ray.dir, d2);
  float8 det = Dot (d1, PVec);

  if (abs (det) > Small_float8)
    {
    float8 inv_det = 1.0 / det;
    VECTOR TVec = Ray.start - p1;
    float8 u = Dot (TVec, PVec) * inv_det;
    if ((u>=0.0) && (u<1.0))
      {
      VECTOR QVec = Cross (TVec, d1);
      float8 v = Dot (Ray.dir, QVec) * inv_det;
      if ((v>=0.0) && ((u+v)<1.0))
        {
        float8 d = Dot (d2, QVec) * inv_det;
        if ((d>Ray.mindist) && (d<dist))
          {
          found = true;
          dist = d;
          ru = u; rv = v;
          }
        }
      }
    }

  d1 = p2 - p4;
  d2 = p3 - p4;
  PVec = Cross (Ray.dir, d2);
  det = Dot (d1, PVec);

  if (abs (det) > Small_float8)
    {
    float8 inv_det = 1.0 / det;
    VECTOR TVec = Ray.start - p4;
    float8 u = Dot (TVec, PVec) * inv_det;
    if ((u>=0.0) && (u<1.0))
      {
      VECTOR QVec = Cross (TVec, d1);
      float8 v = Dot (Ray.dir, QVec) * inv_det;
      if ((v>=0.0) && ((u+v)<1.0))
        {
        float8 d = Dot (d2, QVec) * inv_det;
        if ((d>Ray.mindist) && (d<dist))
          {
          found = true;
          dist = d;
          ru = 1-v; rv = 1-u;
          }
        }
      }
    }

  return found;
  }

} // namespace RAYPP

#endif
