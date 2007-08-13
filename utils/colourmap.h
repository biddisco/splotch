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

#ifndef RAYPP_COLOURMAP_H
#define RAYPP_COLOURMAP_H

#include "kernel/colour.h"
#include "kernel/handle.h"

namespace RAYPP {

class CMAP_ENTRY
  {
  protected:
    float4 minval, maxval;
    float4 fract (float4 value) const
      { return (value-minval)/(maxval-minval); }

  public:
    CMAP_ENTRY () {}
    CMAP_ENTRY (float4 min, float4 max)
      : minval (min), maxval (max) {}
    virtual ~CMAP_ENTRY () {}

    virtual COLOUR Get_Colour (float8) const = 0;

    bool Is_Inside (float8 value) const
      {
      return ((value >= minval) && (value <= maxval));
      }

    void Set_Values (float8 min, float8 max)
      {
      minval = min; maxval = max;
      }
  };

class COLOURMAP
  {
  private:
    vector<HANDLE<CMAP_ENTRY> > Entry;

  public:
    COLOURMAP() {}
    COLOURMAP(COLOUR Col1, COLOUR Col2);

    COLOUR Get_Colour (float8) const;

    void Add_Entry (const HANDLE<CMAP_ENTRY> &);
  };

class UNIFORM_CMAP_ENTRY : public CMAP_ENTRY
  {
  protected:
    COLOUR Colour;

  public:
    UNIFORM_CMAP_ENTRY () {};
    UNIFORM_CMAP_ENTRY (float8 Start, float8 End, const COLOUR& Col)
      : CMAP_ENTRY (Start, End), Colour (Col) {}

    virtual COLOUR Get_Colour (float8) const 
      { return Colour; }
  };

class LINEAR_CMAP_ENTRY: public CMAP_ENTRY
  {
  protected:
    COLOUR Colour1, Colour2;

  public:
    LINEAR_CMAP_ENTRY ()
      : Colour1 (0,0,0), Colour2 (0,0,0) {}
    LINEAR_CMAP_ENTRY (float8 Start, float8 End, 
      const COLOUR &c1, const COLOUR &c2)
      : CMAP_ENTRY (Start, End), Colour1 (c1), Colour2 (c2) {}

    virtual COLOUR Get_Colour (float8 value) const
      {
      return Colour1 + fract(value) * (Colour2 - Colour1);
      }
  };

} // namespace RAYPP

#endif
