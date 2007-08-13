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

#ifndef RAYPP_CAMERA_H
#define RAYPP_CAMERA_H

#include "config/config.h"
#include "kernel/initable.h"
#include "kernel/transformable.h"
#include "kernel/colour.h"

namespace RAYPP {

      /** \class CAMERA kernel/camera.h kernel/camera.h
      Prototype for all camera classes. */
class CAMERA: public INITABLE, public TRANSFORMABLE
  {
  public:
         /*! Returns the light intensity received by the camera from the
         direction described by \a u and \a v. Both \a u and \a v must lie in
         the range [0; 1]. For a normal camera, \a u = \a v = 0 would mean the
         upper left corner and \a u = \a v = 1 the lower right corner of the
         image. */
    virtual COLOUR Calc_Intensity (float8 u, float8 v) const = 0;
  };

} // namespace RAYPP

#endif
