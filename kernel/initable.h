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

#ifndef RAYPP_INITABLE_H
#define RAYPP_INITABLE_H

#include "config/config.h"
#include "kernel/error.h"

/** The namespace enclosing the whole Ray++ library
 */
namespace RAYPP {

/**
   \class INITABLE kernel/initable.h kernel/initable.h
   A class implementing a simple invariant.
 */
class INITABLE
  {
  protected:
    /*!
       This variable indicates if the object is in a valid state (one could
       call it the invariant of the object).

       If <em>false</em>, only functions for the
       configuration of the object may be called (in RAY++, their names usually
       begin with "Set_"). <br>
       If <em>true</em>, you may call only functions
       that leave the state of the object unaltered (that means you can only
       'ask questions' of the object).

       Any constructor must set this variable to <em>false</em>. <br>
       The function Init() must set this variable to <em>true</em>.
      */
    bool initialized;

    /*!
       abbreviation for "Check if Initialized".
       Makes sure that the object is initialized;
       otherwise, emit an error
      */
    void ci () const
      { if (!initialized) error ("Call only allowed after Init()"); }
    /**
       abbreviation for "Check if Not Initialized".
       makes sure that the object has not yet been initialized;
       otherwise, emit an error
      */
    void cni () const
      { if (initialized) error ("Call only allowed before Init()"); }

  public:
    /**
       constructor: sets initialized to false
      */
    INITABLE () : initialized (false) {}
    /** */
    virtual ~INITABLE () {}

    /**
       dummy init function: sets initialized to true.
       This is a placeholder for more useful initialization functions
      */
    virtual void Init () { initialized = true; }
    /**
       dummy deinit function: sets {\tt initialized} to false.
       This is a placeholder for more useful de-initialization functions
      */
    virtual void Deinit () { initialized = false; }
  };

} // namespace RAYPP

#endif
