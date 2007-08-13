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

#ifndef RAYPP_VOLUME_H
#define RAYPP_VOLUME_H

#include "config/config.h"
#include "kernel/initable.h"
#include "kernel/transformable.h"
#include "kernel/colour.h"

namespace RAYPP {

class RAY;
class SHADING_INFO;

/**
  \class VOLUME kernel/volume.h
  Prototype for all volume shaders.
  This class contains functions for calculating the modification of light
  travelling through a medium.
*/
class VOLUME: public INITABLE, public TRANSFORMABLE
  {
  public:
    /*!
      Returns the Index of Refraction for the given Info.
     */
    virtual float8 Refractive_Index (const SHADING_INFO &Info) const = 0;

    /*! */
    virtual COLOUR Calc_new_Importance (const RAY &Ray) const = 0;

    /*! */
    virtual COLOUR Calc_Modified_Colour (const RAY &Ray,
      const COLOUR &Ingoing) const = 0;

    /*! */
    virtual COLOUR Get_Attenuated_Light (const RAY &Ray,
      const COLOUR &Ingoing) const = 0;
  };

} // namespace RAYPP

#endif
