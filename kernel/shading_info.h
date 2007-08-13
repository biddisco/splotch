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

#ifndef RAYPP_SHADING_INFO_H
#define RAYPP_SHADING_INFO_H

#include "config/config.h"
#include "kernel/vector.h"
#include "kernel/colour.h"
#include "kernel/handle.h"
#include "kernel/surface.h"
#include "kernel/inside_info.h"

namespace RAYPP {

/**
  \class SHADING_INFO kernel/shading_info.h kernel/shading_info.h
  Helper class for shading calculations.
*/
class SHADING_INFO
  {
  public:
    /*! */
    TMPHANDLE<SURFACE> Surf;
// evtl. u,v
    /*! */
    VECTOR Intersect_Point;
    /*! */
    VECTOR Incident_Dir;
    /*! */
    VECTOR Normal;
    /*! */
    INSIDE_INFO Ins1, Ins2;
    /*! */
    COLOUR Importance;
    /*! */
    uint1 diffuse_level, specular_level;
  };

} // namespace RAYPP

#endif
