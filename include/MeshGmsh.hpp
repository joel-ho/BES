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
#ifndef MESH_GMSH_HPP
#define MESH_GMSH_HPP

#include <cmath>

#include "../include/Config.hpp"
#include "../include/Mesh.hpp"

class FaceGmshOne: public Face {
  public:
    virtual void initialize();
    virtual void computeCenter ();
    virtual void computeNormal ();
};

class FaceGmshTwo: public Face {
  public:
    virtual void initialize();
    virtual void computeCenter ();
    virtual void computeNormal ();
};

class ElementGmshTwo: public Element {
  public:
    virtual void initialize();
    virtual void computeVolume();
    virtual void computeCentroid();
};

class ElementGmshFour: public Element {
  public:
    virtual void initialize();
    virtual void computeVolume();
    virtual void computeCentroid();
};

#endif