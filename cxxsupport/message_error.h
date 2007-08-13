/*
 * Copyright (c) 2002-2005 Max-Planck-Society
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#ifndef PLANCK_MESSAGE_ERROR_H
#define PLANCK_MESSAGE_ERROR_H

#include <exception>
#include <iostream>
#include <string>

#if defined (PLANCK_STACKTRACE)
#include <execinfo.h>
#endif

inline void show_stackframe()
  {
#if defined (PLANCK_STACKTRACE)
  void *trace[16];
  int trace_size = backtrace(trace, 16);
  char **messages = backtrace_symbols(trace, trace_size);
  std::cerr << "[bt] Execution path:" << std::endl;
  for (int i=0; i<trace_size; ++i)
    std::cerr << "[bt] " << messages[i] << std::endl;
#endif
  }


class Message_error
  {
  private:
    std::string msg;

  public:
    Message_error()
      : msg (std::string("Unspecified error"))
      { std::cerr<<msg<<std::endl; show_stackframe(); }

    explicit Message_error(const std::string &message)
      : msg (message) { std::cerr<<msg<<std::endl; show_stackframe(); }

    virtual const char* what() const
      { return msg.c_str(); }

    virtual ~Message_error() {}
  };

#if defined (PLANCK_CHECKS)

#define PLANCK_DIAGNOSIS_BEGIN try {
#define PLANCK_DIAGNOSIS_END \
} \
catch (Message_error &e) \
  { std::cerr << "Planck exception: " << e.what() << std::endl; throw; } \
catch (std::exception &e) \
  { std::cerr << "std::exception: " << e.what() << std::endl; throw; } \
catch (...) \
  { std::cerr << "Unknown exception" << std::endl; throw; }

#else

#define PLANCK_DIAGNOSIS_BEGIN
#define PLANCK_DIAGNOSIS_END

#endif

#endif
