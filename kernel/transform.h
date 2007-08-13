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

#ifndef RAYPP_TRANSFORM_H
#define RAYPP_TRANSFORM_H

#include "config/config.h"
#include "kernel/vector.h"

namespace RAYPP {

/**
  \class TRANSMAT kernel/transform.h kernel/transform.h
  Helper class for TRANSFORM and STRANSFORM.
*/
class TRANSMAT
  {
  public:
    union
      {
      float4 entry[3][4];
      float4 p[12];
      };

    /*! */
    TRANSMAT () {}
    /*! */
    TRANSMAT (float4, float4, float4,
              float4, float4, float4,
              float4, float4, float4,
              float4, float4, float4);

    /*! */
    TRANSMAT &operator*= (const TRANSMAT &b);
    /*! */
    TRANSMAT &operator+= (const TRANSMAT &b);
    /*! */
    TRANSMAT &operator-= (const TRANSMAT &b);

    /*! */
    TRANSMAT &operator*= (float8 factor);

    /*! */
    TRANSMAT Inverse () const;
    /*! */
    void Invert ();
    /*! */
    void SetToIdentity ();
    /*! */
    void SetToZero ();
    /*! */
    void Transpose ();

    /*! */
    bool Orthogonal () const;
    /*! */
    bool Orthonormal () const;
    /*! */
    bool Scaled_Orthonormal (float8 &factor) const;
    /*! */
    bool Diagonal () const;

    /*! */
    friend ostream &operator<< (ostream &os, const TRANSMAT &mat);
  };

/**
  \class TRANSFORM kernel/transform.h kernel/transform.h
  A class for linear 3D transformations.
*/
class TRANSFORM 
  {
  private:
    TRANSMAT matrix;
    TRANSMAT inverse;

    friend class STRANSFORM;

  public:
    TRANSFORM ()
      {
      matrix.SetToIdentity();
      inverse.SetToIdentity();
      }
    /*! */
    const TRANSMAT &Matrix  () const
      { return matrix; }
    /*! */
    const TRANSMAT &Inverse () const
      { return inverse; }

    /*! */
    void Invert()
      {
      swap(matrix,inverse);
      }
    /*! */
    void Make_Scaling_Transform (const VECTOR &vec);
    /*! */
    void Make_Translation_Transform (const VECTOR &vec);
    /*! */
    void Make_Rotation_Transform (const VECTOR &vec);
    /*! */
    void Make_Axis_Rotation_Transform (const VECTOR &axis, float8 angle);
    /*! */
    void Make_Shearing_Transform
      (float4 xy, float4 xz, float4 yx, float4 yz, float4 zx, float4 zy);

    /*! */
    void Make_General_Transform(const TRANSMAT &trans) 
        { matrix = trans; inverse = matrix.Inverse(); };

    /*! */
    void Add_Transform (const TRANSFORM &trans);  

    /*! */
    bool Orthogonal () const;
    /*! */
    bool Orthonormal () const;
    /*! */
    bool Scaled_Orthonormal (float8 &factor) const;
    /*! */
    bool Diagonal () const;

    /*! */
    VECTOR TransPoint (const VECTOR &vec) const;
    /*! */
    VECTOR InvTransPoint (const VECTOR &vec) const;
    /*! */
    VECTOR TransDirection (const VECTOR &vec) const;
    /*! */
    VECTOR InvTransDirection (const VECTOR &vec) const;
    /*! */
    VECTOR TransNormal (const VECTOR &vec) const;
    /*! */
    VECTOR InvTransNormal (const VECTOR &vec) const;

    /*! */
    friend TRANSFORM Scaling_Transform (const VECTOR &vec);
    /*! */
    friend TRANSFORM Translation_Transform (const VECTOR &vec);
    /*! */
    friend TRANSFORM Rotation_Transform (const VECTOR &vec);
    /*! */
    friend TRANSFORM Axis_Rotation_Transform
      (const VECTOR &axis, float8 angle);
    /*! */
    friend TRANSFORM Shearing_Transform
      (float4 xy, float4 xz, float4 yx, float4 yz, float4 zx, float4 zy);

    /*! */
    friend ostream &operator<< (ostream &os, const TRANSFORM &trans);
  };  

/**
  \class STRANSFORM kernel/transform.h kernel/transform.h
  A compact version of TRANSFORM.
*/
class STRANSFORM 
  {
  private:
    TRANSMAT inverse;

  public:
    /*! */
    STRANSFORM () {inverse.SetToIdentity();}

    /*! */
    STRANSFORM (const TRANSFORM &trans) : inverse (trans.inverse) {}
    /*! */
    operator TRANSFORM () const;

    /*! */
    const TRANSMAT &Inverse () const
      { return inverse; }

    /*! */
    void Add_Transform (const STRANSFORM &strans);

    /*! */
    bool Orthogonal () const;
    /*! */
    bool Orthonormal () const;
    /*! */
    bool Scaled_Orthonormal (float8 &factor) const;
    /*! */
    bool Diagonal () const;

    /*! */
    VECTOR InvTransPoint (const VECTOR &vec) const;
    /*! */
    VECTOR InvTransDirection (const VECTOR &vec) const;
    /*! */
    VECTOR TransNormal (const VECTOR &vec) const;

    /*! */
    friend ostream &operator<< (ostream &, const STRANSFORM &);
  };  

} // namespace RAYPP

#endif
