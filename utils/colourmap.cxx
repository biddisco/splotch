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

#include "utils/colourmap.h"

namespace RAYPP {

COLOURMAP::COLOURMAP(COLOUR Col1, COLOUR Col2)
  {
  Add_Entry(new LINEAR_CMAP_ENTRY (0.0, 1.0, Col1, Col2));
  }

COLOUR COLOURMAP::Get_Colour (float64 val) const
  {
  vector<HANDLE_RAYPP<CMAP_ENTRY> >::const_iterator i;
    for (i=Entry.begin(); i<Entry.end(); ++i)
    {
    if ((*i)->Is_Inside (val))
      return ((*i)->Get_Colour (val));
    }
  cerr << "COLOURMAP: Didn't find colourmap entry" << endl;
  return COLOUR (0,0,0);
  }

void COLOURMAP::Add_Entry (const HANDLE_RAYPP<CMAP_ENTRY> &newentry)
  {
  Entry.push_back (newentry);
  }

} // namespace RAYPP
