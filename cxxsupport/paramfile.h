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

#ifndef PLANCK_PARAMFILE_H
#define PLANCK_PARAMFILE_H

#include <map>
#include <string>
#include <iostream>
#include "cxxutils.h"

class paramfile
  {
  private:
    std::map<std::string,std::string> params;
    bool verbose;

    std::string get_valstr(const std::string &key) const
      {
      std::map<std::string,std::string>::const_iterator loc=params.find(key);
      if (loc!=params.end()) return loc->second;
      throw Message_error ("Error: Cannot find the key \"" + key + "\".");
      }

  public:
    paramfile (const std::string &filename, bool verbose_=true)
      : verbose(verbose_)
      { parse_file (filename, params); }

    bool param_present(const std::string &key) const
      { return (params.find(key)!=params.end()); }

    template<typename T> T find (const std::string &key) const
      {
      T result;
      stringToData(get_valstr(key),result);
      if (verbose)
        std::cout << "Parser: " << key << " = " << dataToString(result)
                  << std::endl;
      return result;
      }
    template<typename T> T find
      (const std::string &key, const T &deflt) const
      {
      if (param_present(key)) return find<T>(key);
      if (verbose)
        std::cout << "Parser: " << key << " = " << dataToString(deflt)
                  << " <default>" << std::endl;
      return deflt;
      }
  };

#endif
