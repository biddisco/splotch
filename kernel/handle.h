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

#ifndef RAYPP_HANDLE_H
#define RAYPP_HANDLE_H

#include "config/config.h"

namespace RAYPP {

/**
  \class HANDLE kernel/handle.h kernel/handle.h
  Pointer with reference counting.
  This class implements a 'smart pointer' with reference count and
  automatic deallocation after deletion of the last reference.
*/
template <typename T> class HANDLE
  {
  private:
    class HANDLE_DATA
      {
      private:
// prevent copying
        HANDLE_DATA (const HANDLE_DATA &);
        HANDLE_DATA &operator=(const HANDLE_DATA &);

      public:
        T *data;
        uint4 counter;

        HANDLE_DATA (T *val)
          : data (val), counter (1) {}
        ~HANDLE_DATA ()
          { delete data; }
      };

    HANDLE_DATA *HandlePtr;

    void dispose()
      { if (HandlePtr) if (--(HandlePtr->counter) == 0) delete HandlePtr; }

  public:
    class BAD_OP {};

    /*!
       Default constructor. Handle points to 0.
      */
    HANDLE () 
      : HandlePtr (0) {}
    /*!
       Handle points to the same object as orig.
      */
    HANDLE (const HANDLE &orig)
      : HandlePtr (orig.HandlePtr)
      { if (HandlePtr) ++(HandlePtr->counter); }
    /*!
       Handle is initialized with contents; reference count is set to 1.
       \warning contents must not point to an automatic variable!
       After calling this constructor, this pointer should not be used again in
       any form! <br>
       Recommended usage of this constructor:
       \code HANDLE<A> foo (new A); \endcode
       Using it this way, one can make no mistakes.
      */
    HANDLE (T *contents)
      {
       if (contents)
         HandlePtr = new HANDLE_DATA (contents);
       else
         HandlePtr = 0;
      }
    /*!
       Destructor. If this is the last reference to the object pointed to,
       it is deleted.
      */
    ~HANDLE ()
      { dispose(); }

    /*!
       The handle is unbound from the object it formerly pointed to
       (this object is eventually destroyed) and bound to the object
       pointed to by orig.
      */
    HANDLE &operator= (const HANDLE &orig)
      {
      if (orig.HandlePtr) ++(orig.HandlePtr->counter);
      dispose();
      HandlePtr = orig.HandlePtr;
      return *this;
      }    
    /*!
       The handle is unbound from the object it formerly pointed to
       (this object is eventually destroyed) and bound to the object
       pointed to by \a contents. <br>
       \warning \a contents must not point to an automatic variable!
       After calling this function, this pointer should not be used again in
       any form! <br>
       Recommended usage of this function:
       \code HANDLE<A> foo; foo = new A; \endcode
       using it this way, one can make no mistakes.
      */
    HANDLE &operator= (T *contents)
      {
      dispose();
      if (contents) HandlePtr = new HANDLE_DATA (contents);
      else HandlePtr = 0;
      return *this;
      }
    /*!
       The handle is unbound from the object it formerly pointed to
       (this object is eventually destroyed) and set to 0.
      */
    void Clear ()
      { dispose(); HandlePtr = 0; }

    /*!
       Returns a pointer to the handle's object.
       If the handle points to 0, program execution is stopped.
      */
    T *operator-> () const
      {
      if (HandlePtr) return HandlePtr->data;
      throw BAD_OP();
      }
    /*!
       Returns a reference to the handle's object.
       If the handle points to 0, program execution is stopped.
      */
    T &operator* () const
      {
      if (HandlePtr) return (*(HandlePtr->data));
      throw BAD_OP();
      }

    /*!
       Returns false, if the handle points to 0, else true.
      */
    bool Valid () const
      { return HandlePtr; }
    /*!
       Returns false, if the handle points to 0, else true.
      */
    operator bool () const
      { return (HandlePtr!=0); }

    T* Pointer() const
      {
      if (HandlePtr) return (HandlePtr->data);
      return 0;
      }

    void swap(HANDLE<T>& other)
      { ::std::swap(HandlePtr,other.HandlePtr); }
  };

namespace std {

template<typename T>
  inline void swap(RAYPP::HANDLE<T>& a, RAYPP::HANDLE<T>& b)
    { a.swap(b); }

}

/**
   \class TMPHANDLE kernel/handle.h kernel/handle.h
   Reference-like handle class.
*/
template <typename T> class TMPHANDLE
  {
  private:
    T *Ptr;

  public:
    /*! */
    TMPHANDLE () 
      : Ptr (0) {}
    /*! */
    TMPHANDLE (const TMPHANDLE &orig)
      : Ptr (orig.Ptr) {}
    /*! */
    TMPHANDLE (const HANDLE<T> orig)
      : Ptr (orig.Pointer()) {}
    /*! */
    ~TMPHANDLE ()
      {}
    /*! */
    TMPHANDLE &operator= (const HANDLE<T> orig)
      { Ptr = orig.Pointer(); return *this; }    
    /*! */
    TMPHANDLE &operator= (const TMPHANDLE orig)
      { Ptr = orig.Ptr; return *this; }    
    /*! */
    TMPHANDLE &operator= (T *contents)
      { Ptr = contents; return *this; }
    /*! */
    T *operator-> () const
      { return Ptr; }
    /*! */
    T &operator* () const
      { return (*Ptr); }
    /*! */
    bool Valid () const
      { return Ptr; }
    /*! */
    operator bool () const
      { return (Ptr!=0); }
  };

} // namespace RAYPP

#endif
