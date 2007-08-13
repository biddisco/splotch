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

#ifdef _MSC_VER
namespace std {
#include <math.h>
}
#else
#include <cmath>
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <queue>
#include <utility>
#include <map>
#include <sstream>

#include "config/typeselect.h" 
 
/** The namespace enclosing the whole Ray++ library 
 */ 
namespace RAYPP {

using namespace ::std;
#ifdef _MSC_VER
template<class T> T    min (const T &x,const T &y) { return (x<y ? x : y);   }
template<class T> T    max (const T &x,const T &y) { return (x>y ? x : y);   }
template<class T> T    abs (const T& v)            { return (v>0 ? v : -v ); }
#if _MSC_VER < 1300 
template<class T> void swap(T &x,T &y) { T tmp = x; x = y; y = tmp; } 
#endif
#endif


/// floating point type (min. 4 bytes, IEEE754 single)
typedef TYPESELECT::GET_FLOAT<4>::RET float4;
/// floating point type (min. 8 bytes, IEEE754 double)
typedef TYPESELECT::GET_FLOAT<8>::RET float8;

/// unsigned byte type (exactly 1 byte)
typedef unsigned char byte;

/// unsigned integer type (min. 1 byte)
typedef TYPESELECT::GET_UINT<1>::RET uint1;
///   signed integer type (min. 1 byte)
typedef TYPESELECT::GET_INT<1>::RET int1;
///   unsigned integer type (min. 2 bytes)
typedef TYPESELECT::GET_UINT<2>::RET uint2;
///   signed integer type (min. 2 bytes)
typedef TYPESELECT::GET_INT<2>::RET int2;
///   unsigned integer type (min. 4 bytes)
typedef TYPESELECT::GET_UINT<4>::RET uint4;
///   signed integer type (min. 4 bytes)
typedef TYPESELECT::GET_INT<4>::RET int4;

class ENDIAN_TEST
  {
  private:
    bool big_end;

  public:
    ENDIAN_TEST()
      {
      const uint2 a=1;
      big_end = (reinterpret_cast<const byte *>(&a))[0] == 0;
      }
    operator bool() const { return big_end; }
  };

const ENDIAN_TEST big_endian;

} // namespace RAYPP

#endif
