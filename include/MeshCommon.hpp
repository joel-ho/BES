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
#ifndef MESH_COMMON_HPP
#define MESH_COMMON_HPP

#include <cmath>
#include <algorithm>
#include <iterator>

#include "../include/FileReader.hpp"
#include "../include/Mesh.hpp"

using namespace std;

class Vertex;
class Face;
class Element;

namespace MeshMath {
  
  extern double computeDistance(
    double *coordsOne, double *coordsTwo, const short &nDim);
  
  extern double computeDistance(
    Vertex *vertexOne, Vertex *vertexTwo, const short &nDim);
  
  extern double computeMagnitude(
    double *vector, const short &nDim);
  
  extern void nonDimVector(double *vector, const short &nDim, double area);
  
  extern void computeTwoDimOutwardNormal (
    double *result, Vertex *vertexOne, Vertex *vertexTwo);
  
  extern double dotProduct (double *vecOne, double *vecTwo, const short &nDim);
  
  extern void crossProduct(
    double *result, Vertex *vertexOne, Vertex *vertexTwo, Vertex *vertexThree);
    
  extern void computeVector(
    Vertex* vertexOne, Vertex* vertexTwo, double* resultantVector, const short &nDim);
  
}

#endif // MESH_COMMON_HPP
