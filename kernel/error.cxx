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

#include <cstring>
#include <cerrno>

#include "kernel/error.h"
#include "kernel/globals.h"

namespace RAYPP {

void message (const string &text)
  {
  Message_Stream << text << endl;
  }

void warning (const string &text)
  {
  Message_Stream << "warning: " << text << endl;
  }

void error (const string &text)
  {
  Error_Stream << "Error: " << text << endl;
  abort();
  }

void fatal (const string &text)
  {
  Error_Stream << "Fatal error: " << text << endl;
  abort();
  }

void internal (const string &text)
  {
  Error_Stream << "Internal error: " << text << endl;
  abort();
  }

void syserror ()
  {
  Error_Stream << "System error: " << strerror (errno) << endl;
  exit (errno);
  }

} // namespace RAYPP
