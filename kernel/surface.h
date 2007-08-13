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

#ifndef RAYPP_SURFACE_H
#define RAYPP_SURFACE_H

#include "config/config.h"
#include "kernel/initable.h"
#include "kernel/transformable.h"
#include "kernel/colour.h"

namespace RAYPP {

class SHADING_INFO;
class FULL_SHADING_INFO;
class INCIDENT_ARRAY;

/**
  \class SURFACE kernel/surface.h kernel/surface.h
  Surface properties.
  This class contains functions for determining the emitted and transmitted
  colour of a surface. 
*/
class SURFACE: public INITABLE, public TRANSFORMABLE
  {
  public:
    /*! */
    virtual void Get_Full_Shading_Info (const SHADING_INFO &In,
      FULL_SHADING_INFO &Out) const = 0;

    /*! */
    virtual COLOUR Get_Total_Importance (const FULL_SHADING_INFO &Info,
      const VECTOR & Dir) const = 0;

    /*! */
    virtual VECTOR Get_MC_Reflected_Dir (const FULL_SHADING_INFO &Info)
      const = 0;

    /*! */
    virtual VECTOR Get_MC_Refracted_Dir (const FULL_SHADING_INFO &Info)
      const = 0;

    /*! */
    virtual VECTOR Get_MC_Diffuse_Dir (const FULL_SHADING_INFO &Info)
      const = 0;

    /*! */
    virtual COLOUR Get_Emitted_Light (const FULL_SHADING_INFO &Info,
      const INCIDENT_ARRAY &Light, const COLOUR &Reflected,
      const COLOUR &Refracted, const INCIDENT_ARRAY &MC_Reflected,
      const INCIDENT_ARRAY &MC_Refracted, const INCIDENT_ARRAY &MC_Diffuse)
      const = 0;

    /*!
      This function must return the colour that is transmitted by the surface
      under the conditions described in Info. Ingoing is the incident colour.
     */
    virtual COLOUR Get_Transmitted_Light 
      (const SHADING_INFO &Info, const COLOUR &Ingoing) const = 0;
  };

} // namespace RAYPP

#endif
