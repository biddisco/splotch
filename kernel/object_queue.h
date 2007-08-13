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

#ifndef RAYPP_OBJECT_QUEUE_H
#define RAYPP_OBJECT_QUEUE_H

#include "config/config.h"
#include "kernel/object.h"
#include "kernel/handle.h"

namespace RAYPP {

/** \class oqentry kernel/object_queue.h kernel/object_queue.h
  Building block for OBJECT_QUEUE.
*/
class oqentry
  {
  public:
    /*! */
    float8 first;
    /*! */
    TMPHANDLE<OBJECT> second;

    oqentry() {}                                    //MR: why?
    /*! */
    oqentry (float8 f, const TMPHANDLE<OBJECT> &s)
      : first(f), second(s) {}

    /*! */
    bool operator< (const oqentry &other) const
// YES, that's correct, since we need the NEAREST intersection first
      { return (first > other.first); }
  };

/** \class OBJECT_QUEUE kernel/object_queue.h kernel/object_queue.h
  Helper class for intersection calculations.
  Inherits priority_queue<oqentry>.
*/
class OBJECT_QUEUE: public priority_queue <oqentry>
  {};

} // namespace RAYPP

#endif
