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

#ifndef RAYPP_GLOBALS_H
#define RAYPP_GLOBALS_H

#include "config/config.h"
#include "kernel/world.h"
#include "kernel/renderer.h"
#include "kernel/handle.h"
#include "kernel/msg_stream.h"

namespace RAYPP {

extern HANDLE<WORLD>    World;    /*!< */
extern HANDLE<RENDERER> Renderer; /*!< */

extern MSG_STREAM Message_Stream; /*!< */
extern MSG_STREAM     Log_Stream; /*!< */
extern MSG_STREAM   Debug_Stream; /*!< */
extern MSG_STREAM   Error_Stream; /*!< */

} // namespace RAYPP

#endif
