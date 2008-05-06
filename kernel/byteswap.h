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

#ifndef RAYPP_BYTESWAP_H
#define RAYPP_BYTESWAP_H

namespace RAYPP {

template<int size> inline void byteswap_helper (byte *)
  {
  bool Error_unspecialized_template_used[-(size*size)];
  // compile time error would be nice here ...
  }
template<> inline void byteswap_helper<1> (byte *)
  {}
template<> inline void byteswap_helper<2> (byte *val)
  {
  swap (val[0],val[1]);
  }
template<> inline void byteswap_helper<4> (byte *val)
  {
  swap (val[0],val[3]); swap (val[1],val[2]);
  }
template<> inline void byteswap_helper<8> (byte *val)
  {
  swap (val[0],val[7]); swap (val[1],val[6]);
  swap (val[2],val[5]); swap (val[3],val[4]);
  }

template<typename T> inline void byteswap (T& val)
  {
  byteswap_helper<sizeof(T)> (reinterpret_cast<byte *> (&val));
  }

} // namespace RAYPP

#endif
