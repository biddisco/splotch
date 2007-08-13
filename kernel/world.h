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

#ifndef RAYPP_WORLD_H
#define RAYPP_WORLD_H

#include "config/config.h"
#include "kernel/initable.h"
#include "kernel/colour.h"

namespace RAYPP {

class RAY;
class VECTOR;
class INSIDE_INFO;
class SHADING_INFO;
class LIGHT_ARRAY;

/**
  \class WORLD kernel/world.h kernel/world.h
  A class holding the scene description.
*/
class WORLD: public INITABLE
  {
  public:
    /*! */
    virtual void Get_Surrounding_Volume
      (const VECTOR &Loc, INSIDE_INFO &result) const = 0;

    /*! */
    virtual bool Get_Next_Intersection
      (const RAY &Ray, float8 &dist, SHADING_INFO &result) const = 0;

    /*! */
    virtual void Get_Lights
      (const VECTOR &Loc, LIGHT_ARRAY &result) const = 0;

    /*! */
    virtual COLOUR Get_Background
      (const VECTOR &Dir) const = 0;
  };

} // namespace RAYPP

#endif
