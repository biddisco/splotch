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

#ifndef RAYPP_RPN_FEXPR_H
#define RAYPP_RPN_FEXPR_H

#include "kernel/kernel.h"
#include "utils/noise.h"

namespace RAYPP {

class RPN_FEXPR
  {
  private:
    enum { maxstack = 100, maxdata=101 };

    typedef enum { FLOAT, FUNCTION } TOKEN_TYPE;

    typedef enum { ADD, SUB, MUL, DIV, SIN, COS, ABS, NEG, V_MAG, V_SCALE,
                   V_ADD, CLAMP, NOISE, V_NOISE, POP, DUP, GET, PUT }
      FUNCTION_TYPE;

    class FUNCMAP: public map<string, FUNCTION_TYPE>
      {
      public:
        FUNCMAP();
      };

    static FUNCMAP funcmap;
    static ::RAYPP::NOISE Noise;

    class EXPR_ENTRY
      {
      public:
        float8 f;
        int4 par;
        FUNCTION_TYPE fct;
        TOKEN_TYPE type;

        EXPR_ENTRY (FUNCTION_TYPE Fct) : fct (Fct), type(FUNCTION) {}
        EXPR_ENTRY (FUNCTION_TYPE Fct, int4 Par)
          : par(Par), fct(Fct), type(FUNCTION) {}
        EXPR_ENTRY (float8 F) : f(F), type (FLOAT) {}
      };

    vector<EXPR_ENTRY> entries;

    float8 data[maxdata];

    void Parse (const string &val);

    inline void stackcheck (int &depth, int before, int after) const;
    inline void datacheck (int pos) const;

    void Validate ();

  public:
    RPN_FEXPR () {}
    RPN_FEXPR (const string &expr)
      { Parse (expr); }

    void Set_Expr (const string &expr)
      { Parse (expr); }

    void Eval ();
    float8 &operator[] (int4 pos)
      { return data[pos]; }
  };

} // namespace RAYPP

#endif
