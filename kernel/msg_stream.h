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

#ifndef RAYPP_MSG_STREAM_H
#define RAYPP_MSG_STREAM_H

#include "config/config.h"

namespace RAYPP {

/**
  \class MSG_STREAM kernel/msg_stream.h
  A utility for redirecting output streams.
*/
class MSG_STREAM
  {
  private:
    ostream *ost;

  public:
    /*! */
    MSG_STREAM ()
      : ost (0) {}

    /*! */
    MSG_STREAM (ostream &s)
      : ost (&s) {}

    /*! */
    void SetTo (ostream &s)
      { ost = &s; }

    /*! */
    void SetToNull ()
      { ost = 0; }

    template<typename T> MSG_STREAM &operator<< (const T &data)
      {
      if (ost) (*ost) << data;
      return *this;
      }
    MSG_STREAM &operator<< (ostream &(data)(ostream&))
      {
      if (ost) (*ost) << data;
      return *this;
      }
  };

} // namespace RAYPP

#endif
