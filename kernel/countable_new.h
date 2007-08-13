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

/* This code is derived from the ideas presented by Kevlin Henney
   in http://www.boost.org/more/count_bdy.htm */

#ifndef RAYPP_COUNTABLE_NEW_H
#define RAYPP_COUNTABLE_NEW_H

#ifdef RAYPP_FAST_HANDLES

#include <cstddef>

namespace RAYPP {

class COUNTABLE_T {};
extern const COUNTABLE_T countable;

class BAD_COUNTABLE_PTR {};

enum { alignment = 4 };
enum { countable_pad = alignment/sizeof(::std::size_t) };

enum { countable_magic_number = 574432 };

inline ::std::size_t &countable_counter (const void *ptr)
  {
  return *(const_cast< ::std::size_t *>
    (static_cast<const ::std::size_t *>(ptr) - countable_pad)); }

inline void countable_check (const void *ptr)
  {
  if (countable_counter(ptr) != countable_magic_number)
    throw BAD_COUNTABLE_PTR();
  }

inline void *operator new (size_t size, const COUNTABLE_T &)
  {
  size_t *allocated = static_cast< ::std::size_t *>
    (::operator new (countable_pad*sizeof(::std::size_t) + size));
  *allocated = countable_magic_number;
  return allocated + countable_pad;
  }

inline void operator delete (void *ptr, const COUNTABLE_T &)
  { ::operator delete (static_cast< ::std::size_t *>(ptr) - countable_pad); }

using ::operator new;
using ::operator delete;
}

#define Cnew new (::RAYPP::countable) 

#else

#define Cnew new

#endif

#endif
