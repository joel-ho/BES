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
#ifndef FACE_HPP
#define FACE_HPP

#include "../include/Options.hpp"
#include "../include/Mesh.hpp"

using namespace std;

class Vertex;
class Face;
class Element;

class Face {
  // Normal vector points out of ownerElement into neighborElement
  
  protected:
    
    short nDim_;
    int nVertex_;
    unsigned long id_;
    double area_;
    double *nonDimNormal_, *dimNormal_;
    
    BoundaryType boundaryType_;
    
    Vertex *center_;
    
    // Pointers to external entities in the mesh 
    // (do not call delete on object pointed to)
    Vertex **vertices_;
    Element *ownerElement_, *neighborElement_;
    
    void initializeBase_(const int &nVertex, const short &nDim);
  
  public:
    
    Face();
    virtual ~Face();
    
    virtual void initialize() = 0;
    
    virtual void computeCenter () = 0;
    virtual void computeNormal () = 0;
    
    void setId(unsigned long id);
    void setSingleVertex(Vertex *vertex, const int &iVertex);
    void setVertices(Vertex **vertices);
    void setBoundaryType(BoundaryType boundaryType);
    void setOwnerElement(Element *ownerElement);
    void setNeighborElement(Element *neighborElement);
    
    void flipNormal();
    
    BoundaryType getBoundaryType();
    
    int queryBoundaryExternal();
    void queryBoundaryExternal(int &boundaryExternal);
    
    Element* getOwnerElement();
    void getOwnerElement(Element *ownerElement);
    
    Element* getNeighborElement();
    void getNeighborElement(Element *neighborElement);
    
    double* getDimNormal();
    // getDimNormal(double *normal) ;
    double getDimNormal(const short &iDim);
    
    double* getNonDimNormal();
    double getNonDimNormal(const short &iDim);
    
    Vertex* getCenter();
    
    double getArea();
    void getArea(double &area);
    
    unsigned long getId();
    
    void printDimNormal();
    void printNonDimNormal();
    
};

#endif // FACE_HPP
