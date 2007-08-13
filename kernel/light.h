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

#ifndef RAYPP_LIGHT_H
#define RAYPP_LIGHT_H

#include "config/config.h"
#include "kernel/initable.h"
#include "kernel/transformable.h"

namespace RAYPP {

class LIGHT_ARRAY;
class VECTOR;

/**
  \class LIGHT kernel/light.h
  Prototype for all light sources.
*/
class LIGHT: public INITABLE, public TRANSFORMABLE
  {
  public:
    /*!
       If this lightsource can illuminate the location Pos, fill the 
       necessary information into Arr.
     */
    virtual void Cast_Light
      (const VECTOR &Pos, LIGHT_ARRAY &Arr) const = 0;
  };

} // namespace RAYPP

#endif
