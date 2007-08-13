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

#ifndef RAYPP_FULL_SHADING_INFO_H
#define RAYPP_FULL_SHADING_INFO_H

#include "config/config.h"
#include "kernel/vector.h"
#include "kernel/colour.h"
#include "kernel/shading_info.h"

namespace RAYPP {

/**
  \class FULL_SHADING_INFO kernel/full_shading_info.h kernel/full_shading_info.h
  An extended SHADING_INFO.
*/
class FULL_SHADING_INFO: public SHADING_INFO
  {
  public:
    /*! */
    bool MC_diffuse : 1;
    /*! */
    bool reflect    : 1;
    /*! */
    bool MC_reflect : 1;
    /*! */
    bool refract    : 1;
    /*! */
    bool MC_refract : 1;

    /*! */
    VECTOR Bumped_Normal;
    /*! */
    VECTOR Reflected;
    /*! */
    VECTOR Refracted;
    /*! */
    COLOUR Reflected_Imp;
    /*! */
    COLOUR Refracted_Imp;
    /*! */
    COLOUR Diffuse_Imp;
  };

} // namespace RAYPP

#endif
