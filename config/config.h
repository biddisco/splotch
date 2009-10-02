/*
 *  Ray++ - Object-oriented ray tracing library
 *  Copyright (C) 1998-2004 Martin Reinecke and others.
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

#ifndef RAYPP_CONFIG_H
#define RAYPP_CONFIG_H

#include <cmath>
#include <vector>
#include <cstdlib>
#ifdef VSS
#include "cxxsupport/datatypes.h"
#else
#include "datatypes.h"
#endif


/** The namespace enclosing the whole Ray++ library 
 */ 
namespace RAYPP {

using namespace ::std;

typedef unsigned char byte;

class ENDIAN_TEST
  {
  private:
    bool big_end;

  public:
    ENDIAN_TEST()
      {
      const uint16 a=1;
      big_end = (reinterpret_cast<const byte *>(&a))[0] == 0;
      }
    operator bool() const { return big_end; }
  };

const ENDIAN_TEST big_endian;

} // namespace RAYPP

#endif
