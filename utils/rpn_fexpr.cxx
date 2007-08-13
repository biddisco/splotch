/*
 *  Ray++ - Object-oriented ray tracing library
 *  Copyright (C) 1998-2004 Martin Reinecke and others.
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

#include "utils/rpn_fexpr.h"

namespace RAYPP {

//temporary
bool get_float(const string &s, float8 &f)
  {
  istringstream tmp (s);
  return (tmp >> f);
  }
bool get_uint(const string &s, int4 &i)
  {
  istringstream tmp (s);
  return (tmp >> i);
  }
bool get_func(const string &s, string &tok, int4 &arg)
  {
  size_t p=0;
  while ((p<s.length()) && (!isdigit(s[p])))
    ++p; 
  tok=s.substr(0,p);
  if (p==s.length()) { arg=1; return true; }
  size_t p2=p;
  while ((p2<s.length()) && (isdigit(s[p2])))
    ++p2;
  istringstream tmp (s.substr(p,p2-p));
  tmp >> arg;
  return (p2==s.length());
  }

// static
RPN_FEXPR::FUNCMAP RPN_FEXPR::funcmap;

// static
NOISE RPN_FEXPR::Noise;

RPN_FEXPR::FUNCMAP::FUNCMAP()
  {
  FUNCMAP &fm = *this;
  fm["+"] = ADD;
  fm["-"] = SUB;
  fm["*"] = MUL;
  fm["/"] = DIV;
  fm["sin"] = SIN;
  fm["cos"] = COS;
  fm["abs"] = ABS;
  fm["neg"] = NEG;
  fm["clamp"] = CLAMP;
  fm["noise"] = NOISE;
  fm["v_mag"] = V_MAG;
  fm["v_scale"] = V_SCALE;
  fm["v_add"] = V_ADD;
  fm["v_noise"] = V_NOISE;
  fm["pop"] = POP;
  fm["dup"] = DUP;
  fm["get"] = GET;
  fm["put"] = PUT;
  }

void RPN_FEXPR::Parse (const string &val)
  {
  entries.clear();
  istringstream input(val);

  while (input)
    {
    bool done = false;

    string tok_str;
    input >> tok_str;
    if (tok_str == "") done = true;

    if (!done)
      {
      float8 f;
      if (get_float(tok_str, f))
        {
        done = true;
        entries.push_back(f);
        }
      }

    if (!done)
      {
      string tok;
      int4 arg;
      if (get_func (tok_str, tok, arg))
        {
        if (funcmap.find(tok) != funcmap.end())
          {
          done = true;
          FUNCTION_TYPE func = funcmap[tok];
          entries.push_back(EXPR_ENTRY(func, arg));
          }
        }
      }

    if (!done)
      {
      Error_Stream << "found something unusable: " << tok_str << endl;
      error ("unexpected token in function string");
      }
    }

  Validate();
  }

inline void RPN_FEXPR::stackcheck (int &depth, int before, int after) const
  {
  if (depth<before) error ("stack underflow");
  depth += after-before;
  if (depth >= maxstack) error ("stack overflow");
  }

inline void RPN_FEXPR::datacheck (int pos) const
  {
  if ((pos<0) || (pos>=maxdata)) error ("data access error");
  }

void RPN_FEXPR::Validate ()
  {
  int4 depth=0;
  for (size_t m=0; m<entries.size(); ++m)
    {
    if (entries[m].type != FUNCTION)
      stackcheck (depth, 0, 1);
    else
      {
      int4 p = entries[m].par;
      switch (entries[m].fct)
        {
        case ADD: case SUB: case MUL: case DIV:
          stackcheck (depth, 2, 1); break;
        case SIN: case COS: case ABS: case NEG:
          stackcheck (depth, 1, 1); break;
        case V_MAG: case CLAMP: case NOISE:
          stackcheck (depth, 3, 1); break;
        case V_SCALE:
          stackcheck (depth, 4, 3); break;
        case V_ADD:
          stackcheck (depth, 6, 3); break;
        case V_NOISE:
          stackcheck (depth, 3, 3); break;
        case POP:
          stackcheck (depth, p, 0); break;
        case DUP:
          stackcheck (depth, p, 2*p); break;
        case GET:
          datacheck (p);
          stackcheck (depth, 0, 1); break;
        case PUT:
          datacheck (p);
          stackcheck (depth, 1, 0); break;
        default:
          error ("unknown function type"); break;
        }
      }
    }
  }

void RPN_FEXPR::Eval ()
  {
  float8 tmp[maxstack];
  float8 *top = &tmp[0];

  for (size_t m=0; m<entries.size(); ++m)
    {
    switch (entries[m].type)
      {
      case FLOAT: *(++top) = entries[m].f; break;
      case FUNCTION :
        switch (entries[m].fct)
          {
          case ADD: top[-1] += *top; --top; break;
          case SUB: top[-1] -= *top; --top; break;
          case MUL: top[-1] *= *top; --top; break;
          case DIV: top[-1] /= *top; --top; break;
          case ABS: *top = abs(*top); break;
          case SIN: *top = sin(*top); break;
          case COS: *top = cos(*top); break;
          case NEG: *top = -(*top); break;
          case V_MAG:
            top[-2] =
              sqrt (top[-2]*top[-2] + top[-1]*top[-1] + *top * *top);
            top-=2;
            break;
          case V_SCALE:
            top[-1] *= *top; top[-2] *= *top; top[-3] *= *top;
            --top;
            break;
          case V_ADD:
            top[-5] += top[-2]; top[-4] += top[-1]; top[-3] += *top;
            top-=3;
            break;
          case V_NOISE:
            {
            VECTOR tmp2 (top[-2], top[-1], *top);
            tmp2 = Noise.DNoise (tmp2);
            top[-2] = tmp2.x; top[-1] = tmp2.y; *top = tmp2.z;
            }
            break;
          case CLAMP: top[-2] = max (top[-1], min (*top, top[-2]));
            top-=2;
            break;
          case NOISE: top[-2] = Noise.Noise (VECTOR (top[-2], top[-1], *top));
            top-=2;
            break;
          case POP: top-=entries[m].par; break;
          case DUP:
            ++top;
            memcpy (top, top-entries[m].par, entries[m].par*sizeof(float8));
            top+=entries[m].par-1;
            break;
          case GET: *(++top) = data[entries[m].par]; break;
          case PUT: data[entries[m].par] = *(top--); break;
          default: error ("internal error"); break;
          }
        break;
      }
    }
  }

} // namespace RAYPP
