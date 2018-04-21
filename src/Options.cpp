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

#include "../include/Options.hpp"

// Map returns default value of 0 if queried with key not in entry
map < string, BoundaryType > BoundaryTypeStringMap = {
  {"INTERNAL", BoundaryType::INTERNAL},
  {"SLIP_WALL", BoundaryType::SLIP_WALL},
  {"FREESTREAM", BoundaryType::FREESTREAM}
};
  
// Map returns default value of 0 if queried with key not in entry
map < BoundaryType, int > ExternalBoundaryMap = {
    {BoundaryType::INTERNAL, -1},
    {BoundaryType::SLIP_WALL, 1},
    {BoundaryType::FREESTREAM, 1},
};