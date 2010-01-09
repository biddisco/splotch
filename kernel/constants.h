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

#ifndef RAYPP_CONSTANTS_H
#define RAYPP_CONSTANTS_H

#include "config/config.h"

const float64 Pi = float64 (3.14159265358979323846264338327950288);
      /*!< is 3.14159265358979323846264338327950288 */

const float32 Small_float32 = float32 (1.0e-4);  /*!< is 1e-4 */
const float32 Large_float32 = float32 (1.0e+4);  /*!< is 1e+4 */
const float32  Huge_float32 = float32 (1.0e+10);  /*!< is 1e+10 */

const float64 Small_float64 = float64 (1.0e-6);  /*!< is 1e-6 */
const float64 Large_float64 = float64 (1.0e+6);  /*!< is 1e+6 */
const float64  Huge_float64 = float64 (1.0e+20);  /*!< is 1e+20 */

const float64 Small_dist   = float64 (1.0e-7);  /*!< is 1e-7 */
const float64 Large_dist   = float64 (1.0e+7);  /*!< is 1e+7 */

#endif
