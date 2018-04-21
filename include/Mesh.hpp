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
#ifndef MESH_HPP
#define MESH_HPP

#include <stdlib.h>
#include <stdexcept>
#include <iostream>
#include <string>
#include <vector>
#include <set>

#include "../include/Options.hpp"
#include "../include/Config.hpp"
#include "../include/MeshCommon.hpp"
#include "../include/Vertex.hpp"
#include "../include/Face.hpp"
#include "../include/Element.hpp"

using namespace std;

class Mesh {
  
  private:
    
    short nDim_, nActiveBoundary_;
    unsigned long nVertex_, nFace_, nElement_;
    
    BoundaryType *activeBoundaryType_;
    
    Vertex **vertices_;
    Face **faces_;
    Element **elements_;
    
    void readGmsh(string meshFilePath);
    
    void saveFaceElement_(vector <Face*> &faces, vector <Element*> &elements);
    void saveActiveBoundary_(set <BoundaryType> &activeBoundaryType);
    void completeInitialization_();
    
    void checkMesh_();
    
  public:
    
    Mesh();
    ~Mesh();
    
    void initialize(Config *config);
    
    short getNumDim();
    short getNumActiveBoundary();
    BoundaryType getActiveBoundaryType(const short &iBoundary);
    
    Vertex* getVertex(const unsigned long &iVertex);
    
    Face* getFace(const unsigned long &iElement);
    
    Element* getElement(const unsigned long &iElement);
    void getElement(Element *resultantElement, const unsigned long &iElement);
    
    unsigned long getNumVertex();
    unsigned long getNumFace();
    unsigned long getNumElement();
    
    void printCoords();
  
};

#endif // MESH_HPP
