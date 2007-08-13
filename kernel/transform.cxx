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

#include "kernel/transform.h"
#include "kernel/constants.h"

namespace RAYPP {

TRANSMAT:: TRANSMAT (float4 a00, float4 a01, float4 a02,
                     float4 a10, float4 a11, float4 a12,
                     float4 a20, float4 a21, float4 a22,
                     float4 a30, float4 a31, float4 a32)
  {
  entry[0][0] = a00;
  entry[1][0] = a01;
  entry[2][0] = a02;

  entry[0][1] = a10;
  entry[1][1] = a11;
  entry[2][1] = a12;

  entry[0][2] = a20;
  entry[1][2] = a21;
  entry[2][2] = a22;

  entry[0][3] = a30;
  entry[1][3] = a31;
  entry[2][3] = a32;
  }

TRANSMAT &TRANSMAT::operator*= (const TRANSMAT &b)
  {
  TRANSMAT a(*this);
  int i,j;
  for (i = 0 ; i < 4 ; ++i)
    for (j = 0 ; j < 3 ; ++j) 
      entry[j][i] = a.entry[0][i] * b.entry[j][0]
                  + a.entry[1][i] * b.entry[j][1]
                  + a.entry[2][i] * b.entry[j][2];
  entry[0][3] += b.entry[0][3];
  entry[1][3] += b.entry[1][3];
  entry[2][3] += b.entry[2][3];
  return *this;
  }  

TRANSMAT &TRANSMAT::operator+= (const TRANSMAT &b)
  {
  for (int i = 0 ; i < 12 ; ++i) p[i] += b.p[i];
  return *this;
  }

TRANSMAT &TRANSMAT::operator-= (const TRANSMAT &b)
  {
  for (int i = 0 ; i < 12 ; ++i) p[i] -= b.p[i];
  return *this;
  }

TRANSMAT &TRANSMAT::operator*= (float8 factor)
  {
  for (int i = 0 ; i < 12 ; ++i) p[i] *= factor;
  return *this;
  }

TRANSMAT TRANSMAT::Inverse () const
  {
  TRANSMAT tmp;
  float8 d;

  tmp.entry[0][0] = entry[1][1]*entry[2][2] - entry[2][1]*entry[1][2];
  tmp.entry[0][1] = entry[0][1]*entry[2][2] - entry[2][1]*entry[0][2];
  tmp.entry[0][2] = entry[0][1]*entry[1][2] - entry[1][1]*entry[0][2];

  tmp.entry[1][0] = entry[1][0]*entry[2][2] - entry[2][0]*entry[1][2];
  tmp.entry[1][1] = entry[0][0]*entry[2][2] - entry[2][0]*entry[0][2];
  tmp.entry[1][2] = entry[0][0]*entry[1][2] - entry[1][0]*entry[0][2];

  tmp.entry[2][0] = entry[1][0]*entry[2][1] - entry[2][0]*entry[1][1];
  tmp.entry[2][1] = entry[0][0]*entry[2][1] - entry[2][0]*entry[0][1];
  tmp.entry[2][2] = entry[0][0]*entry[1][1] - entry[1][0]*entry[0][1];

  d = 1.0 / (entry[0][0]*tmp.entry[0][0] -
             entry[1][0]*tmp.entry[0][1] +
             entry[2][0]*tmp.entry[0][2]);

  if (abs(d) > Huge_float4)
    {
    cerr << "degenerate matrix in TRANSMAT::Inverse()" << endl;
    exit(1);
    }

  tmp.entry[0][0] *= d;
  tmp.entry[2][0] *= d;
  tmp.entry[1][1] *= d;
  tmp.entry[0][2] *= d;
  tmp.entry[2][2] *= d;

  d = -d;

  tmp.entry[1][0] *= d;
  tmp.entry[0][1] *= d;
  tmp.entry[2][1] *= d;
  tmp.entry[1][2] *= d;

  tmp.entry[0][3] = - (tmp.entry[0][0]*entry[0][3] +
                       tmp.entry[0][1]*entry[1][3] +
                       tmp.entry[0][2]*entry[2][3]);
  tmp.entry[1][3] = - (tmp.entry[1][0]*entry[0][3] +
                       tmp.entry[1][1]*entry[1][3] +
                       tmp.entry[1][2]*entry[2][3]);
  tmp.entry[2][3] = - (tmp.entry[2][0]*entry[0][3] +
                       tmp.entry[2][1]*entry[1][3] +
                       tmp.entry[2][2]*entry[2][3]);

  return tmp;
  }

void TRANSMAT::Invert ()
  {
  *this = Inverse();
  }

void TRANSMAT::SetToIdentity ()
  {
  for (int m=0; m<12; ++m) p[m]=0;
  entry[0][0]=entry[1][1]=entry[2][2]=1;
  }

void TRANSMAT::SetToZero ()
  {
  for (int m=0; m<12; ++m) p[m]=0;
  }

void TRANSMAT::Transpose ()
  {
  swap(entry[0][1], entry[1][0]);
  swap(entry[0][2], entry[2][0]);
  swap(entry[1][2], entry[2][1]);
  entry[0][3] = entry [1][3] = entry[2][3] = 0.0;
  }

bool TRANSMAT::Orthogonal () const
  {
  if (abs (entry[0][0]*entry[1][0] +
           entry[0][1]*entry[1][1] +
           entry[0][2]*entry[1][2]) > Small_float4) return false;
  if (abs (entry[0][0]*entry[2][0] +
           entry[0][1]*entry[2][1] +
           entry[0][2]*entry[2][2]) > Small_float4) return false;
  if (abs (entry[2][0]*entry[1][0] +
           entry[2][1]*entry[1][1] +
           entry[2][2]*entry[1][2]) > Small_float4) return false;
  return true;
  }

bool TRANSMAT::Orthonormal () const
  {
  if (!Orthogonal()) return false;
  if (abs (entry[0][0]*entry[0][0] +
           entry[0][1]*entry[0][1] +
           entry[0][2]*entry[0][2] - 1.0) > Small_float4) return false;
  if (abs (entry[1][0]*entry[1][0] +
           entry[1][1]*entry[1][1] +
           entry[1][2]*entry[1][2] - 1.0) > Small_float4) return false;
  if (abs (entry[2][0]*entry[2][0] +
           entry[2][1]*entry[2][1] +
           entry[2][2]*entry[2][2] - 1.0) > Small_float4) return false;
  return true;
  }

bool TRANSMAT::Scaled_Orthonormal (float8 &factor) const
  {
  if (!Orthogonal()) return false;
  factor = entry[0][0]*entry[0][0] +
           entry[0][1]*entry[0][1] +
           entry[0][2]*entry[0][2];
  if (abs (entry[1][0]*entry[1][0] +
           entry[1][1]*entry[1][1] +
           entry[1][2]*entry[1][2] - factor) > Small_float4)
    return false;
  if (abs (entry[2][0]*entry[2][0] +
           entry[2][1]*entry[2][1] +
           entry[2][2]*entry[2][2] - factor) > Small_float4)
    return false;
  factor = sqrt (factor);
  return true;
  }

bool TRANSMAT::Diagonal () const
  {
  if ((abs (entry[1][0]) > Small_float4) ||
      (abs (entry[2][0]) > Small_float4) ||
      (abs (entry[2][1]) > Small_float4) ||
      (abs (entry[0][1]) > Small_float4) ||
      (abs (entry[0][2]) > Small_float4) ||
      (abs (entry[1][2]) > Small_float4)) return false;

  return true;
  }

ostream &operator<< (ostream &os, const TRANSMAT &mat)
  {
  for (int i=0;i<4;++i)
    os << mat.entry[0][i] << ' ' << mat.entry[1][i]
       << ' ' << mat.entry[2][i] << endl;
  return os;
  }

void TRANSFORM::Make_Scaling_Transform (const VECTOR &vec)
  {
  if ((vec.x<Small_float4) || (vec.y<Small_float4) || (vec.z<Small_float4))
    {
    cerr << "TRANSFORM: invalid scaling transformation" << endl;
    exit(1);
    }

  matrix.SetToIdentity();
  matrix.entry[0][0]=vec.x;
  matrix.entry[1][1]=vec.y;
  matrix.entry[2][2]=vec.z;

  inverse.SetToIdentity();
  inverse.entry[0][0]=1.0/vec.x;
  inverse.entry[1][1]=1.0/vec.y;
  inverse.entry[2][2]=1.0/vec.z;
  }

void TRANSFORM::Make_Translation_Transform (const VECTOR &vec)
  {
  matrix.SetToIdentity();
  matrix.entry[0][3]=vec.x;
  matrix.entry[1][3]=vec.y;
  matrix.entry[2][3]=vec.z;

  inverse.SetToIdentity();
  inverse.entry[0][3]=-vec.x;
  inverse.entry[1][3]=-vec.y;
  inverse.entry[2][3]=-vec.z;
  }

void TRANSFORM::Make_Rotation_Transform (const VECTOR &vec)
  {
  TRANSMAT tmp;
  VECTOR Radian_Vector = vec*Pi/180.0;
  float8 cosx, cosy, cosz, sinx, siny, sinz;

  matrix.SetToIdentity();
  cosx = cos (Radian_Vector.x);
  sinx = sin (Radian_Vector.x);
  cosy = cos (Radian_Vector.y);
  siny = sin (Radian_Vector.y);
  cosz = cos (Radian_Vector.z);
  sinz = sin (Radian_Vector.z);

  matrix.entry[1][1] =  cosx;
  matrix.entry[2][2] =  cosx;
  matrix.entry[2][1] =  sinx;
  matrix.entry[1][2] = -sinx;
  inverse = matrix;
  inverse.Transpose();

  tmp.SetToIdentity();
  tmp.entry[0][0] =  cosy;
  tmp.entry[2][2] =  cosy;
  tmp.entry[2][0] = -siny;
  tmp.entry[0][2] =  siny;
  matrix *= tmp;
  tmp.Transpose();
  inverse *= tmp;

  tmp.SetToIdentity();
  tmp.entry[0][0] =  cosz;
  tmp.entry[1][1] =  cosz;
  tmp.entry[1][0] =  sinz;
  tmp.entry[0][1] = -sinz;
  matrix *= tmp;
  tmp.Transpose();
  inverse *= tmp;
  }

void TRANSFORM::Make_Axis_Rotation_Transform
  (const VECTOR &axis, float8 angle)
  {
  VECTOR V = axis.Norm();
  angle *= Pi/180.0;
  float8 cosx = cos (angle), sinx = sin (angle);

  matrix.SetToZero();
  matrix.entry[0][0] = V.x * V.x + cosx * (1.0 - V.x * V.x);
  matrix.entry[1][0] = V.x * V.y * (1.0 - cosx) + V.z * sinx;
  matrix.entry[2][0] = V.x * V.z * (1.0 - cosx) - V.y * sinx;
  matrix.entry[0][1] = V.x * V.y * (1.0 - cosx) - V.z * sinx;
  matrix.entry[1][1] = V.y * V.y + cosx * (1.0 - V.y * V.y);
  matrix.entry[2][1] = V.y * V.z * (1.0 - cosx) + V.x * sinx;
  matrix.entry[0][2] = V.x * V.z * (1.0 - cosx) + V.y * sinx;
  matrix.entry[1][2] = V.y * V.z * (1.0 - cosx) - V.x * sinx;
  matrix.entry[2][2] = V.z * V.z + cosx * (1.0 - V.z * V.z);
  inverse = matrix;
  inverse.Transpose();
  }

void TRANSFORM::Make_Shearing_Transform
  (float4 xy, float4 xz, float4 yx, float4 yz, float4 zx, float4 zy)
  {
  matrix.SetToIdentity();

  matrix.entry[1][0] = xy;
  matrix.entry[2][0] = xz;
  matrix.entry[0][1] = yx;
  matrix.entry[2][1] = yz;
  matrix.entry[0][2] = zx;
  matrix.entry[1][2] = zy;

  inverse = matrix.Inverse();
  }

void TRANSFORM::Add_Transform (const TRANSFORM &trans)
  {
  matrix *= trans.matrix;
  TRANSMAT tmp=trans.inverse;
  tmp*=inverse;
  inverse = tmp;
  }

bool TRANSFORM::Orthogonal () const
  {
  return matrix.Orthogonal();
  }

bool TRANSFORM::Orthonormal () const
  {
  return matrix.Orthonormal();
  }

bool TRANSFORM::Scaled_Orthonormal (float8 &factor) const
  {
  return matrix.Scaled_Orthonormal (factor);
  }

bool TRANSFORM::Diagonal () const
  {
  return matrix.Diagonal();
  }

VECTOR TRANSFORM::TransPoint (const VECTOR &vec) const
  {
  const float4 *p = matrix.p;
  return VECTOR
    (vec.x*p[0] + vec.y*p[1] + vec.z*p[2] + p[3],
     vec.x*p[4] + vec.y*p[5] + vec.z*p[6] + p[7],
     vec.x*p[8] + vec.y*p[9] + vec.z*p[10]+ p[11]);
   }

VECTOR TRANSFORM::InvTransPoint (const VECTOR &vec) const
  {
  const float4 *p = inverse.p;
  return VECTOR
    (vec.x*p[0] + vec.y*p[1] + vec.z*p[2] + p[3],
     vec.x*p[4] + vec.y*p[5] + vec.z*p[6] + p[7],
     vec.x*p[8] + vec.y*p[9] + vec.z*p[10]+ p[11]);
  }

VECTOR TRANSFORM::TransDirection (const VECTOR &vec) const
  {
  const float4 *p = matrix.p;
  return VECTOR
    (vec.x*p[0] + vec.y*p[1] + vec.z*p[2],
     vec.x*p[4] + vec.y*p[5] + vec.z*p[6],
     vec.x*p[8] + vec.y*p[9] + vec.z*p[10]);
  }

VECTOR TRANSFORM::InvTransDirection (const VECTOR &vec) const
  {
  const float4 *p = inverse.p;
  return VECTOR
    (vec.x*p[0] + vec.y*p[1] + vec.z*p[2],
     vec.x*p[4] + vec.y*p[5] + vec.z*p[6],
     vec.x*p[8] + vec.y*p[9] + vec.z*p[10]);
  }

VECTOR TRANSFORM::TransNormal (const VECTOR &vec) const
  {
  const float4 *p = inverse.p;
  return VECTOR
    (vec.x*p[0] + vec.y*p[4] + vec.z*p[8],
     vec.x*p[1] + vec.y*p[5] + vec.z*p[9],
     vec.x*p[2] + vec.y*p[6] + vec.z*p[10]);
  }

VECTOR TRANSFORM::InvTransNormal (const VECTOR &vec) const
  {
  const float4 *p = matrix.p;
  return VECTOR
    (vec.x*p[0] + vec.y*p[4] + vec.z*p[8],
     vec.x*p[1] + vec.y*p[5] + vec.z*p[9],
     vec.x*p[2] + vec.y*p[6] + vec.z*p[10]);
  }

TRANSFORM Scaling_Transform (const VECTOR &vec)
  {
  TRANSFORM trans;
  
  trans.Make_Scaling_Transform (vec);
  return trans;
  }

TRANSFORM Translation_Transform (const VECTOR &vec)
  {
  TRANSFORM trans;
  
  trans.Make_Translation_Transform (vec);
  return trans;
  }

TRANSFORM Rotation_Transform (const VECTOR &vec)
  {
  TRANSFORM trans;
  
  trans.Make_Rotation_Transform (vec);
  return trans;
  }

TRANSFORM Axis_Rotation_Transform (const VECTOR &axis, float8 angle)
  {
  TRANSFORM trans;
  
  trans.Make_Axis_Rotation_Transform (axis, angle);
  return trans;
  }

TRANSFORM Shearing_Transform
  (float4 xy, float4 xz, float4 yx, float4 yz, float4 zx, float4 zy)
  {
  TRANSFORM trans;
  
  trans.Make_Shearing_Transform (xy, xz, yx, yz, zx, zy);
  return trans;
  }

ostream &operator<< (ostream &os, const TRANSFORM &t)
  {
  os << "Transform\n{\n" << t.matrix << t.inverse << "}" << endl;
  return os;
  }


STRANSFORM::operator TRANSFORM () const
  {
  TRANSFORM trans;
  trans.matrix = inverse.Inverse();
  trans.inverse = inverse;
  return trans;
  }

void STRANSFORM::Add_Transform (const STRANSFORM &strans)
  {
  TRANSMAT tmp=strans.inverse;
  tmp*=inverse;
  inverse = tmp;
  }

bool STRANSFORM::Orthogonal () const
  {
  return inverse.Orthogonal();
  }

bool STRANSFORM::Orthonormal () const
  {
  return inverse.Orthonormal();
  }

bool STRANSFORM::Scaled_Orthonormal (float8 &factor) const
  {
  if (!inverse.Scaled_Orthonormal (factor)) return false;
  factor = 1.0/factor;
  return true;
  }

bool STRANSFORM::Diagonal () const
  {
  return inverse.Diagonal();
  }

VECTOR STRANSFORM::InvTransPoint (const VECTOR &vec) const
  {
  const float4 *p = inverse.p;
  return VECTOR (vec.x*p[0] + vec.y*p[1] + vec.z*p[2] + p[3],
                 vec.x*p[4] + vec.y*p[5] + vec.z*p[6] + p[7],
                 vec.x*p[8] + vec.y*p[9] + vec.z*p[10]+ p[11]);
  }

VECTOR STRANSFORM::InvTransDirection (const VECTOR &vec) const
  {
  const float4 *p = inverse.p;
  return VECTOR (vec.x*p[0] + vec.y*p[1] + vec.z*p[2],
                 vec.x*p[4] + vec.y*p[5] + vec.z*p[6],
                 vec.x*p[8] + vec.y*p[9] + vec.z*p[10]);
  }

VECTOR STRANSFORM::TransNormal (const VECTOR &vec) const
  {
  const float4 *p = inverse.p;
  return VECTOR (vec.x*p[0] + vec.y*p[4] + vec.z*p[8],
                 vec.x*p[1] + vec.y*p[5] + vec.z*p[9],
                 vec.x*p[2] + vec.y*p[6] + vec.z*p[10]);
  }

ostream &operator<< (ostream &os, const STRANSFORM &t)
  {
  os << "STransform\n{\n" << t.inverse << "}" << endl;
  return os;
  }

} // namespace RAYPP
