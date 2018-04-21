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
#ifndef VERTEX_HPP
#define VERTEX_HPP

#include "../include/Options.hpp"
#include "../include/MeshCommon.hpp"

using namespace std;

class Vertex {
  
  private:
    
    unsigned long id_;
    double *coords_;
    
  public:
    
    Vertex();
    ~Vertex();
    
    void initialize(const short &nDim);
    void initialize(double* coords, const short &nDim);
    
    void setId(const unsigned long &id);
    
    void setCoords(double *coords, const short &nDim);
    void setSingleCoord(double singleCoord, const short &iCoord);
    
    unsigned long getId();
    
    void getCoords(double* resultantCoords, const short &nDim);
    void getSingleCoord(double &singleCoords, const short &iCoord);
    double getSingleCoord(const short &iCoord);
    
    void printCoords(const short &nDim);
  
};

#endif // VERTEX_HPP
