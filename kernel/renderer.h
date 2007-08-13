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

#ifndef RAYPP_RENDERER_H
#define RAYPP_RENDERER_H

#include "config/config.h"
#include "kernel/initable.h"
#include "kernel/colour.h"
#include "kernel/light_array.h"

namespace RAYPP {

class RAY;

/**
  \class RENDERER kernel/renderer.h kernel/renderer.h
  Prototype for all renderers.
*/
class RENDERER: public INITABLE
  {
  public:
    /*! */
    virtual COLOUR Trace_Ray (const RAY &Ray) const = 0;

    /*! */
    virtual COLOUR Get_Pixel
      (float8 u, float8 v, float8 du, float8 dv) const = 0;
    /*! */
    virtual COLOUR Trace_Camera_Ray (float8 u, float8 v) const = 0;
    /*! */
    virtual void Calc_Illumination (const VECTOR &Loc, const COLOUR &Imp,
      INCIDENT_ARRAY &Arr) const = 0;
  };

} // namespace RAYPP

#endif
