/*
  Basic Euler Solver (BES) - A lightweight Euler Solver for fluid dynamics.
  Copyright (C) 2018, Joel Ho Mun Onn

  BES is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  BES is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with BES.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once
#ifndef FILE_READER_HPP
#define FILE_READER_HPP

// #include <cmath>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <iterator>

using namespace std;

namespace ReadFile {
  
  extern void trimWhiteSpace(string &stringIn);
  
  extern void splitString(
    const string &stringIn, const string &delimiter, string *vectorOut);
  
  extern void splitStringToNum(
    const string &stringIn, const string &delimiter, double *vectorOut);
  
  extern void encodeFaceIdentifier (string &identifier, vector<unsigned long> vertexNumbers);
  
};

#endif // FILE_READER_HPP
