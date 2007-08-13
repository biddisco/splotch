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

#ifndef RAYPP_LIGHT_ARRAY_H
#define RAYPP_LIGHT_ARRAY_H

#include "config/config.h"
#include "kernel/vector.h"
#include "kernel/colour.h"
#include "kernel/inside_info.h"

namespace RAYPP {

/** \class LIGHT_ENTRY kernel/light_array.h kernel/light_array.h
  Building block for LIGHT_ARRAY.
*/
class LIGHT_ENTRY
  {
  public:
    /*! */
    COLOUR Intensity;
    /*! */
    VECTOR Position;
    /*! */
    INSIDE_INFO Inside;

    /*! */
    LIGHT_ENTRY () {}
    /*! */
    LIGHT_ENTRY 
      (const COLOUR &intens, const VECTOR &posit, const INSIDE_INFO &Ins)
      : Intensity (intens), Position (posit), Inside (Ins) {}
  };

/** \class LIGHT_ARRAY kernel/light_array.h kernel/light_array.h
  Helper class for illumination calculations.
  Inherits vector<LIGHT_ENTRY>.
*/
class LIGHT_ARRAY: public vector<LIGHT_ENTRY>
  {
  public:
    /*! */
    COLOUR Ambient;

    /*! */
    LIGHT_ARRAY () : Ambient (0,0,0) {}
  };

/** \class INCIDENT_ENTRY kernel/light_array.h kernel/light_array.h
  Building block for INCIDENT_ARRAY.
*/
class INCIDENT_ENTRY
  {
  public:
    /*! */
    COLOUR Intensity;
    /*! */
    VECTOR Direction;

    /*! */
    INCIDENT_ENTRY () {}
    /*! */
    INCIDENT_ENTRY 
      (const COLOUR &intens, const VECTOR &direct)
      : Intensity (intens), Direction (direct) {}
  };

/** \class INCIDENT_ARRAY kernel/light_array.h kernel/light_array.h
  Helper class for illumination calculations.
  Inherits vector<INCIDENT_ENTRY>.
*/
class INCIDENT_ARRAY: public vector<INCIDENT_ENTRY>
  {
  public:
    /*! */
    COLOUR Ambient;

    /*! */
    INCIDENT_ARRAY () : Ambient (0,0,0) {}
  };

} // namespace RAYPP

#endif
